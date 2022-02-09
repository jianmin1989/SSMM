/**********************************************************************/
/*                HLS-VIIRS shape fusion                             **/
/*                                                                   **/
/*                                                                   **/                                                       
/*  Version 1                                                        **/
/*  Developer:                                                       **/
/*           Jianmin Wang       South Dakota State University        **/
/*           Xiaoyang Zhang     South Dakota State University        **/
/*                                                                   **/
/*                                                                   **/
/*  Zhang et al., 2020,                                              **/
/*   Development and Evaluation of a NewAlgorithm for Detecting 30m  **/
/*      Land Surface Phenology from VIIRS and HLS Time Series        **/
/*  ISPRS Journal of Photogrammetry and Remote Sensing               **/
/*                                                                   **/   
/**********************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
#include <regex.h>
#include "cproj.h"
#include "proj.h"


#define Row_block 600
#define Col_block 600

void correlation_J1(int *tndvi,int *viirs,float *correl);
void fusionMain_Anci(int *tndvi, long *pixelTile);
void pixel_reproj(double *letftop, int utmzone, long *pixelTile);
char *replaceWord(const char *s, const char *oldW, const char *newW);


short ancievi[Row_block][Col_block][2000]; 
long Row_start_all=-1, Row_end_all=-1, Col_start_all=-1, Col_end_all=-1; 
int SEARCHSIZE=5, Row_start=-1, Row_end=-1, Col_start=-1, Col_end=-1, Nitems=0, nrowblock, ncolblock;
long int Off_set=0;
int i9,j9;
int PERIOD=244, Row=3660, Col=3660, Row_anci=2400, Col_anci=2400, FILL=32767;

main(int argc, char *argv[]){
	//Files
	gzFile fpin[2000],fpanci[2000], *outfile;
	FILE *fpinfo;

	//File names or other strings used
	char fnex_in[1000], fnex_ref[1000], fout[1000], dates[PERIOD][500], tmpname2[500], *tmpname3, *filename;
	char *ret, tmpname[5000];

	//Map info for reprojection
	double letftop[4];
	int utmzone;
	long pixelTile[10], pixelTile1[10], nh=36, nv=18;
	int rows, rowe, cols, cole, irow, jcol;


	//Other information
	int nancitiles=1;
	int tileh[100],tilev[100], tilehnow=-1, tilevnow=-1, ci,cj;
	int tndvi[1000];
	int k, m, n, i, j, kk, buf, k1, k2;
	short h[10], fsndvi[2000];
	
	//Read the paramters
	if((argc<1)||(fpinfo=fopen(argv[1],"r"))==NULL)
	{
		printf("Please prepare an input file following the input_example.txt\n");
		exit(1);
	}

	fscanf(fpinfo,"%s",fout);

	fscanf(fpinfo,"%s",fnex_in);
	fscanf(fpinfo,"%d",&FILL);
	fscanf(fpinfo,"%d",&PERIOD); 
	fscanf(fpinfo,"%lf",&letftop[0]);	//lefttop is for the target imagery
	fscanf(fpinfo,"%lf",&letftop[1]);
	fscanf(fpinfo,"%d",&utmzone); 

	fscanf(fpinfo,"%s",fnex_ref);


	fscanf(fpinfo,"%d",&nancitiles);
	for(i=0;i<nancitiles;i++) fscanf(fpinfo,"%d",&tileh[i]);	//tileh is for the ancillary imagery
	for(i=0;i<nancitiles;i++) fscanf(fpinfo,"%d",&tilev[i]);	//tileh is for the ancillary imagery	
	fscanf(fpinfo,"%d",&SEARCHSIZE);
	for(i=0;i<PERIOD;i++){ 
		fscanf(fpinfo,"%s",tmpname2); 	
		strcpy(dates[i], tmpname2); 
	}	
	fclose(fpinfo);
	//Read the paramters
	//printf("dates[-1]=%s\n", dates[PERIOD-1]);


	if((Row_block <= SEARCHSIZE*2+1)||(Col_block <= SEARCHSIZE*2+1)||(Row_anci <= SEARCHSIZE*2+1)){
		printf("\nPlease use a smaller HALF WINDOW!!!\n");
		exit(1);
	} 




	//Open the output file to write
	if((outfile=gzopen(fout, "wb"))==NULL){
		printf("cannot open the fused time series\n");
		exit(1);
	}


	//Open the input vegetation index
	for(k=0;k<PERIOD;k++){
		filename=replaceWord(fnex_in, "DDDDDDD", dates[k]);
		if((fpin[k]=gzopen(filename,"rb"))==NULL){
			printf("Cann't open file %s\n", filename);
			exit(1);
		}
	}



	//Deal with the data
	buf=2;
	for(i9=0;i9<Row;i9++){	//i9<Row
		for(j9=0;j9<Col;j9++){
			if((i9%10==0)&&(j9==0)&&(i9>0)) printf("line %d over %d lines finished\n", i9, Row);
			for(k=0;k<PERIOD;k++){
				gzread(fpin[k],h, buf);
				tndvi[k]=(int)h[0];
			}

			//Get the map location in VIIRS projection
			pixelTile[0]=i9;    //top to downm-center of HLS
			pixelTile[1]=j9;    //left to right --center of HLS
			for(i=2;i<10;i++) pixelTile[i]=-10000;
			pixel_reproj(letftop, utmzone, pixelTile);   ///12/2018 xyz 
			
			if((pixelTile[4]<0)||(pixelTile[4]>nh*Col_anci-1)||(pixelTile[5]<0)||(pixelTile[5]>nv*Row_anci-1)) 
				goto Skip_fusion;

			//upperleft corners of the required block SEARCHSIZE. 
			pixelTile[6]=pixelTile[4]-SEARCHSIZE;
			pixelTile[7]=pixelTile[5]-SEARCHSIZE;
			if(pixelTile[6]<0) pixelTile[6]=0;
			if(pixelTile[7]<0) pixelTile[7]=0;
			//lowerright corners of the required block SEARCHSIZE. 
			pixelTile[8]=pixelTile[4]+SEARCHSIZE;
			pixelTile[9]=pixelTile[5]+SEARCHSIZE;
			if(pixelTile[8]>(long)60*Col_anci-1) pixelTile[8]=(long)60*Col_anci-1;
			if(pixelTile[9]>(long)20*Row_anci-1) pixelTile[9]=(long)20*Row_anci-1;


			//Read the anci data and do the fusion.
			if((pixelTile[6]>=Col_start_all)&&(pixelTile[7]>=Row_start_all)&&
			(pixelTile[8]<=Col_end_all)&&(pixelTile[9]<=Row_end_all)){
				fusionMain_Anci(tndvi,pixelTile);			
			} 
			else {
				Row_start_all=pixelTile[7];
				Row_end_all=Row_start_all+Row_block-1;
				if(Row_end_all>(long)20*Row_anci-1) Row_end_all=(long)20*Row_anci-1;
				Col_start_all=pixelTile[6];
				Col_end_all=Col_start_all+Col_block-1;
				if(Col_end_all>(long)60*Col_anci-1) Row_end_all=(long)60*Col_anci-1;

						
				rows=0;	//which range to read in ancievi
				pixelTile1[5]=pixelTile[7];
				while(pixelTile1[5]<=Row_end_all){
					cols=0;
					pixelTile1[4]=pixelTile[6];

					while(pixelTile1[4]<=Col_end_all){						
						pixelTile1[0]=pixelTile1[4]/Col_anci;
						pixelTile1[1]=pixelTile1[5]/Row_anci;
						pixelTile1[2]=pixelTile1[4] % Col_anci;
						pixelTile1[3]=pixelTile1[5] % Row_anci;		
				
						if((pixelTile1[0]!=tilehnow)||(pixelTile1[1]!=tilevnow)){
							if((tilehnow>=0)&&(tilevnow>=0))
								for(k=0;k<PERIOD;k++)  gzclose(fpanci[k]);
							Row_start=-1; Row_end=-1; Col_start=-1; Col_end=-1; tilehnow=-1;tilevnow=-1;
							for(kk=0;kk<nancitiles;kk++) {
								if((pixelTile1[0]!=tileh[kk])||(pixelTile1[1]!=tilev[kk])) continue;  
								tilehnow=pixelTile1[0];
								tilevnow=pixelTile1[1];

								sprintf(tmpname2, "%02d", tilehnow); 
								tmpname3=replaceWord(fnex_ref, "HH", tmpname2);
								sprintf(tmpname2, "%02d", tilevnow);
								tmpname3=replaceWord(tmpname3, "VV", tmpname2);
								for(k=0;k<PERIOD;k++){
									filename=replaceWord(tmpname3, "DDDDDDD", dates[k]);
									if((fpanci[k]=gzopen(filename,"rb"))==NULL){
										printf("Cann't open file %s ", filename);
										exit(1);
									} 
									
								}
								break;
							}
							//The range to be read no matter current tile exists or not. 
							Row_start=pixelTile1[3];
							nrowblock=Row_anci-Row_start;
							if(rows+nrowblock>Row_block)
								nrowblock=Row_block-rows;
							Row_end = Row_start+nrowblock-1;
							Col_start=pixelTile1[2];
							ncolblock=Col_anci-Col_start;
							if(cols+ncolblock>Col_block)
								ncolblock=Col_block-cols;
							Col_end = Col_start+ncolblock-1;		
											

						}
						
						if((tilehnow>=0)&&(tilevnow>=0)){
							for(k1=rows;k1<rows+nrowblock;k1++){
								i=Row_start+k1-rows;
								Off_set=((long) 1) *(Col_anci*i+Col_start)*sizeof(short);
								for(k=0;k<PERIOD;k++) gzseek(fpanci[k], Off_set, SEEK_SET);

								for(k2=cols;k2<cols+ncolblock;k2++){
									for(k=0;k<PERIOD;k++){
										gzread(fpanci[k],h, buf);  
										ancievi[k1][k2][k]=h[0];
									}
								
								}
							} 
						
						} else {
							for(k1=rows;k1<rows+nrowblock;k1++){
								for(k2=cols;k2<cols+ncolblock;k2++){
									for(k=0;k<PERIOD;k++){
										ancievi[k1][k2][k]=-FILL;  
									}
								}
							}
						}
						
						cols=cols+ncolblock;
						pixelTile1[4]=pixelTile1[4]+ncolblock;
					} //end while

					rows=rows+nrowblock;
					pixelTile1[5]=pixelTile1[5]+nrowblock;
				} //end while

				
				fusionMain_Anci(tndvi, pixelTile);
				
					
			} //end read new ancievi

			Skip_fusion:
			//write down the results
			for(k=0;k<PERIOD;k++)
				fsndvi[k]=(short)tndvi[k];
			gzwrite(outfile,fsndvi,sizeof(short)*PERIOD);
			//fwrite(fsndvi, sizeof(short), PERIOD, outfile);



		}
	}
		
	//Close files
	gzclose(outfile)	;	
	for(k=0;k<PERIOD;k++)
		gzclose(fpin[k]);
	if((tilehnow>=0)&&(tilevnow>=0))
		for(k=0;k<PERIOD;k++)  gzclose(fpanci[k]);


}



void pixel_reproj(double *letftop, int utmzone, long *pixelTile)
{
long i,j,k,m,n,irow,icol;
 int CMLEN=256;
/* GCTP projection parameters */
char *ptr;
double incoor[2];
double outcoor[2];
long insys = 0;
long inzone;
double inparm[15];
long inunit;
long indatum;
long ipr = 0;
long jpr = 999;
long jprinv = 3;
long outsys;
long outzone;
double outparm[15];
long outunit;
long outdatum;
long proj;
long zonec;

long flg;		/* error flag for conversion of C version   */
long inflag;            /* flag of tile which is in the input image */
char DoAll;             /* flag of selection tile map               */ 

char file27[CMLEN],file83[CMLEN],libgctp[CMLEN],efile[CMLEN],file1[CMLEN];
char inname[CMLEN], outname[CMLEN],headername[CMLEN];
FILE *in,*inptr,*hdrptr,*outptr;

double in_xSize,in_ySize;                         /*Input image pixel size */
double out_xSize,out_ySize;
double in_left,in_top,in_right,in_bottom;         /*Input image range      */
long in_row,in_col;     /* number of row and column of input map      */
long nband;             /* number of band of input map                */
long dsize;             /* Datatype size char=1 short=2 int=4 float=8 */
long offset;            /* offset of input map pointer                */
unsigned char *data;                /* store BIP map data                         */
unsigned char *BackValue;           /* store background data (outside the image)  */

double out_left,out_top,out_right,out_bottom;      /*Image range of each tile         */
double out_corner[4][2];                           /*Four corners of ISG image        */ 
double MinX,MaxX,MinY,MaxY;          /*range of tiles input map covered */
long nCol,nRow;
double inc_x,inc_y;
int tilecolrow[5];

double all_left,all_top,all_right,all_bottom; 
long row,col; 

/* Assign file names */
strcpy (efile,"error_file.txt");
strcpy (file27,"");
strcpy (file83,"");
strcpy (libgctp,"");

// in_left=699960.000000;   ///east
 // in_top=4800000.000000;   //North
// inzone=18;  ///UTMzone input


//input parameters for the input projection
insys=1;  ///UTM
inzone=utmzone;
for(k=0;k<15;k++) 
  inparm[k]=0;
in_xSize=30;
in_ySize=30;
in_col=Col;
in_row=Row;
inunit=2;

in_left=letftop[0];
in_top=letftop[1];


//input parameters for the output projection
nRow=Row_anci;
nCol=Col_anci;
outsys=16;  ///Sinusoidal
outzone=62;
for(k=0;k<15;k++) 
  outparm[k]=0;
out_xSize=463.31271652;
out_ySize=463.31271652;
outunit=2;

out_top=10007670.505136;
out_left=-20015225.181744;

//  Set values to change the datum and radius 
indatum=0;
if ((inparm[0] != 0) && (inparm[0] != 6370997))
   if ((insys != 1) && (insys != 2))
      indatum = -1;

outdatum=0;
if ((outparm[0] != 0) && (outparm[0] != 6370997))
   if ((outsys != 1) && (outsys != 2))
      outdatum = -1;

// for State Plane projection, get data files 
if (insys == 2||outsys == 2) {
   ptr = (char *)getenv("LIBGCTP");
   strcpy(libgctp,ptr);
   sprintf(file27,"%s/nad27sp", libgctp);
   sprintf(file83,"%s/nad83sp", libgctp);
}


m=pixelTile[0];
n=pixelTile[1];
incoor[0]=in_left+n*in_xSize+0.5*in_xSize; 
incoor[1]=in_top-m*in_ySize-0.5*in_ySize;        

// printf("tt2 %f  %f  \n",incoor[0],incoor[1]);

gctp(incoor,&insys,&inzone,inparm,&inunit,&indatum,&ipr,efile,&jpr,file1,
  outcoor,&outsys,&outzone,outparm,&outunit,&outdatum,
  file27,file83,&flg);
//printf("tt3 %f  %f  \n",outcoor[0],outcoor[1]);

// if((m==10)&&(n==10))
// printf("x==%f;y==%f\n", outcoor[0],outcoor[1]);
for(i=0;i<5;i++)  
	tilecolrow[i]=0;



all_left=(outcoor[0]-out_left)/(out_xSize);
all_top=(out_top-outcoor[1])/(out_ySize);
icol=(long)(all_left/nCol);  //TileH
irow=(long)(all_top/nRow);  //TileV
col=(long)((all_left -icol*nCol));
row=(long)((all_top - irow*nRow)); 
if((col>=nCol)||(row>=nRow)) {
	printf("i9=%d, j9=%d, x=%f y=%f all_left=%f all_top=%f tileh=%d tilev=%d col=%d row=%d\n", i9, j9, outcoor[0], outcoor[1], all_left, all_top, icol, irow, col, row);
	exit(1);
}
//	printf("roi1==%d, col==%d\n", row, col);
tilecolrow[0]=icol;  //tile H
tilecolrow[1]=irow;  //Tile V 
tilecolrow[2]=col; //col
tilecolrow[3]=row; //row
//	outcoor[0]=(double)col;
//	outcoor[1]=(double)row;


for(i=0;i<4;i++) pixelTile[i]=tilecolrow[i]; 
pixelTile[4]=(long)all_left;
pixelTile[5]=(long)all_top;
}





// C program to search and replace all occurrences of a word with other word. Scrached from 
//https://www.geeksforgeeks.org/c-program-replace-word-text-another-given-word/ 
// Function to replace a string with another string 
char *replaceWord(const char *s, const char *oldW, const char *newW) { 
	char *result; 
	int i, cnt = 0; 
	int newWlen = strlen(newW); 
	int oldWlen = strlen(oldW); 

	// Counting the number of times old word 
	// occur in the string 
	for (i = 0; s[i] != '\0'; i++) { 
		if (strstr(&s[i], oldW) == &s[i]) { 
			cnt++; 
			i += oldWlen - 1; // Jumping to index after the old word. 
		} 
	} 

	// Making new string of enough length 
	result = (char *)malloc(i + cnt * (newWlen - oldWlen) + 1); 

	i = 0; 
	while (*s) { 
		// compare the substring with the result 
		if (strstr(s, oldW) == s) { 
			strcpy(&result[i], newW); 
			i += newWlen; 
			s += oldWlen; 
		} else
			result[i++] = *s++; 
	} 

	result[i] = '\0'; 
	return result; 
} 





void correlation_J1(int *tndvi,int *viirs,float *correl){

	int i,j,i1,j1,k,lag, laglamda, h;
	float x1,x2,x3,x4, x5, lamda, dlamda, baselam = 1.0;
	int m,n,n1;
	float mean,vimean;
	float x[PERIOD], y[PERIOD];
	float xn[PERIOD], yn[PERIOD];
	float a,c, tx,ty,txy,txx,tyy, meanx, meany,slop,res;
	float singg,r, rr;
	float dxy, dxy1,dxyest, dxymeanx, dyyest,wac,mses,mseu,promeses,rmse, mbe;
	float sumx,sumy,ygmb,ygma, xgmb,xgma,ypred,xpred,rmpds,rmpdu,acs,acu,spdu;
	float spod;
	float max,min;
	short Nfiles=PERIOD;
	short df;
	float sigfi;
	float t;

	//lag=10;
	lag=7;
	//Jianmin 06/2019: added the lamda
	laglamda = 0;
	dlamda = 0.03;
	
	
	//printf("correl[5]=%f ", correl[5]);
	for(h=-laglamda;h<=laglamda;h++){
		lamda=baselam + h*dlamda;
		for(k=-lag;k<lag;k++){
			// m=0;
			n=0;
			for(j=0;j<Nfiles;j++){
				x[j]=FILL;
				y[j]=FILL;   ///01/2019

				j1=round(lamda*(j+k));
				if((j1>0)&&(j1<Nfiles)){ 
					if((viirs[j1]>0)&&(viirs[j1]<FILL)&&(tndvi[j]>0)&&(tndvi[j]<FILL)){
						x[n]=(float)viirs[j1]/10000.0;
						y[n]=(float)tndvi[j]/10000.0;
						n++;
					}
				}
			}



			//******regression
			r=0;
			rr=0;
			rmse=0;
		
			sigfi=0;
			slop=0;
			t=0;
			a=0;

			//if(n>3)
			//if(n>PERIOD/10)
			if(n>PERIOD/20)   /// 01/2019
			{
				tx=0;
				ty=0;
				txy=0;
				txx=0;
				tyy=0;

				for(i=0;i<n;i++)
				{
					tx=tx+x[i];                //*total x **x1*
					ty=ty+y[i];                 //**total y**x2*
					txy=txy+x[i]*y[i];            //**total x*y   **x3
					txx=txx+x[i]*x[i];             //**total x*x   **x4
					tyy=tyy+y[i]*y[i];            //***total y*y   ***x5
				}
				meanx=tx/n;
				meany=ty/n;

				x1=tx;
				x2=ty;
				x3=txy;
				x4=txx;
				x5=tyy;

				/***OLQ****01/2019****
				if((n*x4-x1*x1)==0)
				slop=0;
				else 
				slop=(float)(n*x3-x1*x2)/(n*x4-x1*x1);
				a=(ty-slop*tx)/n;  ///intercept
				*****end OLSQ***/

				singg=(n*txx-tx*tx)*(n*tyy-ty*ty);
				if(singg==0)
					r=0;
				else
					r=(n*txy-tx*ty)/sqrt(singg); //**correlation**


				//****GMFR**01/2019***
				x1=0.0;
				x2=0.0;

				for(i=0;i<n;i++)
				{
					x1=x1+(x[i]-meanx)*(x[i]-meanx);
					x2=x2+(y[i]-meany)*(y[i]-meany);
				}
				if(x1>0)  slop=sqrt(x2/x1);
				if(r<0)  slop=-slop;
				a=meany-slop*meanx;

				//***end GMFR***

			// }

				//******end regression**
				//calculting fitted RMS 01/2019
				ty=0;
				tyy=0;

				for(i=0;i<n;i++){
					ty=a+slop*x[i];
					tyy=tyy+(y[i]-ty)*(y[i]-ty);
				}
				tyy=tyy/n;


				// if(correl[0]<r){  //01/2019
				if((correl[4]>tyy)&&(correl[0]<r)){   
					//printf("lamda=%f k=%d ", lamda, k);
					correl[0]=r;
					correl[1]=a*10000.0;
					correl[2]=slop;
					correl[3]=k;
					correl[4]=tyy;
					correl[5]=lamda;
				}

			}  ///n>PERIOD/10

		} ///lag
	}
	
	
	

} //end




void fusionMain_Anci(int *tndvi,long *pixelTile){   ///05/2019 Modified by Jianmin from XYZ fusionVIIRS_HLS
	int i,j,i1,j1,k,k1=0, k2,anci[PERIOD]; /// 01/2019
	int col,row,col1,row1;
	//int c1=Col_anci, r1=Row_anci;
	int size, c;
	float correl[10],r,a,b,rmse, lamda;  /// 01/2019
	int xstart=0, xend=PERIOD, ystart=0, yend=PERIOD;
	int nvalanci;

	size=SEARCHSIZE;

	col=pixelTile[0]*Col_anci+pixelTile[2]-Col_start_all;
	row=pixelTile[1]*Row_anci+pixelTile[3]-Row_start_all;
	xstart = col-size;
	if(xstart<0) xstart=0;
	ystart = row-size;
	if(ystart<0) ystart=0;
	xend = col+size;
	if(xend > Col_block-1) xend=Col_block-1;
	yend = row+size;
	if(yend > Row_block-1) yend=Row_block-1;


	//printf("\n\ni9=%d, j9=%d, col=%d, row=%d, xstart=%d, xend=%d, ystart=%d, yend=%d, ianci=%d, janci=%d, tilev=%d, tileh=%d\n", i9, j9, col, row, xstart, xend, ystart, yend, (row+Row_start_all)%Col_anci, (col+Col_start_all)%Col_anci, (row+Row_start_all)/Row_anci, (col+Col_start_all)/Col_anci);
	//printf("EVI before fusion\n"); for(k=0;k<PERIOD;k++) printf("%d ", tndvi[k]);printf("\n");
	//printf("Anci evi same location\n");for(k=0;k<PERIOD;k++) printf("%d ", ancievi[row][col][k]);printf("\n");

	i1=-1;
	j1=-1;
	r=0.0;
	a=0.0;
	b=0.0;
	c=0;
	lamda=1.0;
	
	rmse=30000.00;		/// 01/2019
	for(k=0;k<10;k++)
		correl[k]=0.0;
	correl[5]=lamda;	//Jianmin 06/2019
	correl[4]=rmse;		//Jianmin 06/2019 move from fusion function inlizating RMSE selecting the smallest one 01/2019 
	
	
	
	//printf("BBB row=%d Row_start=%d, nrowblock=%d\n", row, Row_start, nrowblock);
	//add this paragraph to say if current pixel works fine, then no window searching
	if(ancievi[row][col][0]==-FILL){
		//printf("Can't file Tile h%02dv%02d for ancillary images!\nFor better prediction of the fused time series, please make it available!\n", (col+Col_start_all)/Col_anci, (row+Row_start_all)/Row_anci);
		goto WINDOWSEARCH1;
	}
	nvalanci=0;
	for(k=0;k<PERIOD;k++){
		//printf("%d ", row-Row_start);
		anci[k]=ancievi[row][col][k];
		if(anci[k] < 10000) nvalanci++;
		
	}
	//printf("nvalanci=%d\n", nvalanci);
	//for(k=0;k<PERIOD;k++) printf("%d ", anci[k]);printf("\n");
	if(nvalanci<PERIOD*0.4) goto WINDOWSEARCH1;
	correlation_J1(tndvi,anci,correl);
	if((r<correl[0])&&(rmse>correl[4])){    
		r=correl[0];  //r
		a=correl[1];  //a
		b=correl[2];  //b
		c=(int)correl[3];  //t (-10 ~10)
		rmse=correl[4];  //rmse  01/2019 
		lamda=correl[5];

		j1=col;
		i1=row;
		k2=k1;
	}	
	//printf("central pixel r=%f\n", r);
	if(r>0.8) goto WINDOWSEARCH2;

	WINDOWSEARCH1:




	
	//printf("Can't file Tile h%02dv%02d for ancillary images!\nFor better prediction of the fused time series, please make it available!\n", (j+Col_start_all)/Col_anci, (i+Row_start_all)/Row_anci);
	
	for(i=ystart;i<=yend;i++){
		for(j=xstart;j<=xend;j++){
			if((i==row)||(j==col)) continue;
			if(ancievi[i][j][0]==-FILL) {
				//printf("\n\nCan't file Tile h%dv%2d for ancillary images!\nFor better prediction of the fused time series, please make it available!\n", (j+Col_start_all)/Col_anci, (i+Row_start_all)/Row_anci);
				continue;
			}
			nvalanci=0;
	  		for(k=0;k<PERIOD;k++){
				anci[k]=ancievi[i][j][k];
				if(anci[k] < 10000) nvalanci++;
			}
			//printf("nvalanci=%d ", nvalanci);
			if(nvalanci<PERIOD*0.4) continue;	//07/2019
			correlation_J1(tndvi,anci,correl);
			if((r<correl[0])&&(rmse>correl[4])){    /// 01/2019
				r=correl[0];  //r
				a=correl[1];  //a
				b=correl[2];  //b
				c=(int)correl[3];  //t (-10 ~10)
				rmse=correl[4];  //rmse  01/2019 
				lamda=correl[5];

				j1=j;
				i1=i;
				k2=k1;
			}
			k1++;
		}
	
	}

	WINDOWSEARCH2:
	//printf("r=%f, a=%f, b=%f, beta=%d, lamda=%f, rmse=%f, kk=%d, rppixeli=%d, rppixelj=%d, i1=%d, j1=%d\n",r,a,b,c,lamda, rmse, k2, pixelTile[3], pixelTile[2], i1, j1);	
	if(r>0.6){
		//if(i9==0) printf(" r=%f", r);
		for(k=0;k<PERIOD;k++){
			//ancievi_print[k]=ancievi[i1][j1][k];  //Jianmin 06/2019: added for print check.
			//ancievi_prediction[k]=FILL; 
			k1=round(lamda*(k+c));
			if((k1>0)&&(k1<PERIOD)){
				//if((i9==1308)&&(j9==81)&&(k==208)) printf("ancievi=%d, k1=%d, a=%f, b=%f, res=%d, k=%d, resfloat=%f, resint=%d\n", ancievi[i1][j1][k1], k1, a, b, ancievi_prediction[k], k, a+b*ancievi[i1][j1][k1], (int)(a+b*ancievi[i1][j1][k1]));
				if((tndvi[k]==FILL)&&(ancievi[i1][j1][k1]!=FILL)&&(ancievi[i1][j1][k1]>0)){ ///01/2019 xyz
				//if((snow[k]>0)&&(ancievi[i1][j1][k1]!=FILL)&&(ancievi[i1][j1][k1]>0)){ ///01/2019 Jianmin   win0non0: win0 SEARCHSIZE=0, non0 SNOW[K]>0
				  	tndvi[k]=a+b*ancievi[i1][j1][k1];
					
				}
			}
			
		}
		
	}
//printf("EVI after fusion\n"); for(k=0;k<PERIOD;k++) printf("%d ", tndvi[k]);printf("\n");
//	
}









