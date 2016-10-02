#include <stddef.h>
#include <stdio.h>
#include <inttypes.h>
#define ARRAYSIZE(x)  (sizeof(x)/sizeof(*(x)))
#define ELEMENT_COUNT(X) (sizeof(X) / sizeof((X)[0]))
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#define pi 3.1416
#define ROWS 240
#define COLS 320
float array[ROWS][COLS], array_fn[ROWS][COLS];    /*Array for building image values*/
float maxradSCH[ROWS][COLS];
float maxblurSCH[ROWS][COLS];
float *xpolar, *ypolar;   /*Array for building polar co-ordinates*/
float *R, *L, *GF, *D, *resmat;			  /*Array for building rmin:rmax values and the gaussian filter*/
float *FINE_SCH;
int *X_th, *Y_th, *X_fn, *Y_fn;

/*FUNCTION 0: Function for creating delay*/
void delay(){
puts("\nEnter any number and press enter\n");
int del=scanf("%d",&del);
}
/*FUNCTION 1: Function to import CSV elements*/

void arraybuild(void) {
   const char filename[] = "sample_image_csv_run.csv";
   FILE *file = fopen(filename, "r"); /*Opens the file*/
   if ( file )
   {
      size_t i, j, k;
      char buffer[4096], *ptr;
      for ( i = 0; fgets(buffer, sizeof buffer, file); ++i ) /*Reads each line from file*/
      {
         for ( j = 0, ptr = buffer; j < ARRAYSIZE(*array); ++j, ++ptr ) /*Parse CSV elements into array*/
         {
            array[i][j] = strtof(ptr, &ptr);
			array_fn[i][j]=array[i][j];
         }
      }
      fclose(file);
	  /* Code for printing (DISABLED!! */
      /*for ( j = 0; j < i; ++j )
      {
         puts("\n");
		 printf("\narray[%lu]: ", (long unsigned)j);
		 puts("\n");
         for ( k = 0; k < ARRAYSIZE(*array); ++k )
         {
            printf("%f ", array[j][k]);
         }
		} */
      putchar('\n');
	}
}

/*FUNCTION 2: Function to dynamically build polar co-ordinates array*/

int plrarrayalloc(float thetorg, float thetain, float thetafin, int elemc) {
	while (thetain<=thetafin)
		{
        thetain=thetain+thetorg;
		elemc=elemc+1;
		}
    xpolar=(float*)malloc(elemc*sizeof(float)); /*Array dynamic memory allocation XPOLAR*/
    ypolar=(float*)malloc(elemc*sizeof(float)); /*Array dynamic memory allocation YPOLAR*/
    /*printf("\nthetain=%f thetafin=%f lasttheta=%f arrsize=%d\n",thetain,thetafin,testz,elemc);*/
return elemc;
}

/*FUNCTION 3: Function to store polar co-ordinates array*/

int plrarraybld(float thetorg, float thetain, float thetafin, int elemc, float Cx, float Cy, float rad) {
int flag=0;
while (thetain<=thetafin)
	{
        float xsub=rad*sin(thetain); float ysub=rad*cos(thetain);
        *(xpolar+elemc)=Cx-xsub; /*printf("%f\n",*(xpolar+elemc));*/
        *(ypolar+elemc)=Cy-ysub; /*printf("%f\n",*(ypolar+elemc));*/

        if (*(xpolar+elemc)>=ROWS || *(xpolar+elemc)<=1 ||
            *(ypolar+elemc)>=COLS || *(ypolar+elemc)<=1)
        {
            flag=1;
        }
        thetain=thetain+thetorg;
        elemc=elemc+1;
}
return flag;
}

/*FUNCTION 3: Function to perform line integral search for PUPIL boundary*/

float pupilintgl(int arraysize) {
int i=0;
float linval=0;
float div=(float)arraysize;
        for (i=0;i<arraysize;i++)
        {
            int xcor=(int)(*(xpolar+i));int ycor=(int)(*(ypolar+i));
            float val=array[xcor][ycor];
            linval=linval+val;
            /*printf("\nxcor=%d, ycor=%d i=%d, val=%f linval=%f\n",xcor,ycor,i,val,linval);*/
        }
float Lfin=(linval/div);
return Lfin;
}

/*FUNCTION 4: Function to perform line integral search for IRIS boundary*/

float irisintgl(int arraysize) {
int i=0;
float linval=0;
float div=(float)arraysize;
int r8=(int)(arraysize/8);
int r3i8=(int)((3*arraysize)/8);
int r5i8=(int)((5*arraysize)/8);
int r7i8=(int)((7*arraysize)/8);

		for(i=1;i<=r8;i++) {
            int xcor=(int)(*(xpolar+i));int ycor=(int)(*(ypolar+i));
            float val=array[xcor][ycor];
            linval=linval+val;
		}

		for(i=(1+r3i8);i<=r5i8;i++) {
            int xcor=(int)(*(xpolar+i));int ycor=(int)(*(ypolar+i));
            float val=array[xcor][ycor];
            linval=linval+val;
		}

		for(i=(1+r7i8);i<arraysize;i++) {
            int xcor=(int)(*(xpolar+i));int ycor=(int)(*(ypolar+i));
            float val=array[xcor][ycor];
            linval=linval+val;
		}
float Lfin=(2*linval)/div;
return Lfin;
}

/*FUNCTION 5: Function to dynamically build rmin:rmax array used in Partiald*/

int rmimxPD(float rmiPD, float rmaPD) {
int i=0;
	while (rmiPD<=rmaPD)
		{
        rmiPD=rmiPD+1;
		i=i+1;
		}
    R=(float*)malloc(i*sizeof(float)); /*Array dynamic memory allocation R*/
	L=(float*)malloc(i*sizeof(float));
return i;
}

/*FUNCTION 6: Function to make rmin:rmax elements used in Partiald*/

void rmimxbdPD(float rmiPD, float rmaPD, int arsz) {
int i=0;
for (i=0;i<arsz;i++)
{
*(R+i)=rmiPD;
rmiPD=rmiPD+1;
}
}

/*FUNCTION 6: Function to make 1-D Gaussian filter*/

void gaussian(int numelems, float sigmaGD){
GF=(float*)malloc(numelems*sizeof(float));
float denm=2*sigmaGD*sigmaGD;
float inttemp=(float)numelems;
float startnum=-((inttemp-1)/2);
int i; float j;
for (i=0;i<numelems;i++)
{
*(GF+i)=startnum+(float)i;
}
/*-----------------------------------------------------------------------------*/
for (i=0;i<numelems;i++)
{
j=*(GF+i);
*(GF+i)=(exp(-((j*j)+(startnum*startnum))/denm))/(pi*denm);
}
/*-----------------------------------------------------------------------------*/
j=0;
for (i=0;i<numelems;i++)
{
j=j+(*(GF+i));
}
/*-----------------------------------------------------------------------------*/
for (i=0;i<numelems;i++)
{
(*(GF+i))=(*(GF+i))/j;
}
/*-----------------------------------------------------------------------------*/
}

/*FUNCTION 7: Function for convolution*/

int convn(int nMat1, int nMat2, float **mat1, float **mat2, float **resmat)
{
  /*mat1=(float*)malloc((nMat1)*sizeof(float));
  mat2=(float*)malloc((nMat2)*sizeof(float));*/
  *resmat=(float*)malloc((nMat1+nMat2-1)*sizeof(float));
  /*--------------------------------------------------------*/
  /*printf("\nEnter the values in Mat1\n");*/
  int iteri;
  /*for (iteri=0; iteri<nMat1; iteri++) {
	scanf("%f",mat1+iteri);
	printf("mat1(%d)=%f\t",iteri,*(mat1+iteri));}
  printf("\nEnter the values in Mat2\n");
  for (iteri=0; iteri<nMat2; iteri++) {
	scanf("%f",mat2+iteri);
	printf("mat1(%d)=%f\t",iteri,*(mat2+iteri));}*/
  /*--------------------------------------------------------*/

  int n=0;
  int SignalLen=nMat1;
  int KernelLen=nMat2;
  /*printf("\n\nValues of ResMAT are\n");*/
  for (n = 0; n < SignalLen + KernelLen - 1; n++)
  {
    int kmin, kmax, k;

    ((*resmat)[n]) = 0;

    kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
    kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

    for (k = kmin; k <= kmax; k++)
    {
      ((*resmat)[n]) += ((*mat1)[k]) * ((*mat2)[n-k]);
    }
    /*printf("\nResmat(%d)=%f",n,((*resmat)[n]));*/
  }
  /*puts("\n");*/

  /*--------------------------------------------------------*/
  /* EQUALISING THE ARRAY DIMENSIONS TO MAT1 */
  int diff=abs((nMat1+nMat2-1)-nMat1);/*printf("\ndiff=%d",diff);*/

  int sublast=diff/2;/*printf("sublast=%d",sublast);*/
  int substart=diff-sublast;/*printf("substart=%d",substart);*/
  int totred=substart+sublast;/*printf("totred=%d",totred);*/
  float *temparr;
  temparr=(float*)malloc(nMat1*sizeof(float));/*printf("mallocsize=%d\n",nMat1);*/
  int counter=0;
  for(iteri=(substart);iteri<((nMat1+nMat2-1)-sublast);iteri++) {
        *(temparr+counter)=((*resmat)[iteri]);
          counter=counter+1;
  }

free(*resmat);
*resmat=(float*)malloc((counter)*sizeof(float));
for(iteri=0;iteri<counter;iteri++) {
        ((*resmat)[iteri])=*(temparr+iteri);
        /*printf("\nResmatfinal(%d)=%f",iteri,*(temparr+iteri));*/
  }
free(temparr);
return(counter);
}

/*FUNCTION 8: Function for finding maximum in an 1D array*/

int arymax1d(float **arrayinp, int array_sz) {

float max_num;
int iteri, fin_idx=1;
max_num=((*arrayinp)[0]);
for (iteri=0; iteri<array_sz; iteri++)
{
	if (((*arrayinp)[iteri])>max_num)
	{
		max_num=((*arrayinp)[iteri]);
		fin_idx=1+iteri;
	}
}
return fin_idx;
}

/*FUNCTION 9: Function for finding maximum in a 2D array*/

float max2d(float a[][320])
{
   int i,j;
   float mymax = a[0][0];
   *(FINE_SCH+0)=0;
   *(FINE_SCH+1)=0;
   for(i = 0; i < ROWS; i++)
      {for(j = 0; j < COLS; j++)
         if(a[i][j] > mymax)
         {
            mymax = a[i][j];
            *(FINE_SCH+0)=i;
            *(FINE_SCH+1)=j;
         }}
	float radret=maxradSCH[(int)*(FINE_SCH+0)][(int)*(FINE_SCH+1)];
    return radret;
}

/*FUNCTION 10: Function for finding minimum in a 2D array*/

/*float min2d(float a[][320])
{
   int i,j;
   float mymin = a[0][0];

   for(i = 0; i < ROWS; i++)
      for(j = 0; j < COLS; j++)
         if(a[i][j] < mymin)
         {
            mymin = a[i][j];
         }
   return mymin;
}
*/

/* LINE INTEGRAL FUNCTION*/

float lineint(float nLI, float rLI, int CxLI, int CyLI, int respLI, int resiLI) {
float Lxx=0;
float theta=(2*pi)/nLI; float thetainval=theta, thetafinval=2*pi;
int flagchk=0; int elmcnt=0;		/*Storing the flag for out of bounds cases, elmcnt for polar array*/
int arrsize=plrarrayalloc(theta,thetainval,thetafinval,elmcnt);  /*Storing the size of polar array*/
float CXLI=(float)CxLI; float CYLI=(float)CyLI;
flagchk=plrarraybld(theta,thetainval,thetafinval,elmcnt,CXLI,CYLI,rLI);
if (flagchk==1){Lxx=0;/*printf("\nFlag returned was\t%d\n",flagchk);*/} /*For any out of bounds-coordinates, no line integral possible*/
else {
        if (respLI==0) { /*Perform line integral search for pupil boundary*/
        Lxx=pupilintgl(arrsize);
        }
        if (resiLI==0) {  /*Perform line integral search for iris boundary*/
        Lxx=irisintgl(arrsize);
        }
} 				   /*Integral code finishes here*/
return Lxx;
}

/* PARTIALD FUNCTION*/

int partiald(int CxPD, int CyPD, float rminPD, float rmaxPD, float sigmaPD, float nPD, int respPD, int resiPD) {
int arrszPD=rmimxPD(rminPD, rmaxPD);
rmimxbdPD(rminPD, rmaxPD, arrszPD);
int i=0;
for (i=0;i<arrszPD;i++)
{
*(L+i)=lineint(nPD, *(R+i), CxPD, CyPD, respPD, resiPD);
if (*(L+i)==0){
*(L+i)=-291;
/*printf("\nBreak Encountered: Flag==1\n");*/
break;
}
}
int j;
int numelmts=i;
/*-------------------------------------------------------------*/
D=(float*)malloc((numelmts)*sizeof(float));
float temp=0;
for (j=0;j<numelmts;j++)
{
*(D+j)=*(L+j)-temp;
/*printf("\n%f\t%f\t%f\tjval=%d\n",*(L+j),temp,*(R+j),j);*/
temp=*(L+j);
/*printf("\n%f at jval=%d\n",*(D+j),j);*/
}
*(D+0)=0;
/*puts("\nDebug 1\n");*/
/*Differernce array *D is built*/
gaussian(5, sigmaPD);
/*Gaussian filter *GF is built*/
int resfinsize=convn(numelmts,5,&D,&GF,&resmat);
int max_indx=arymax1d(&resmat,resfinsize);
max_indx=max_indx-1;
return max_indx;
}

/* SEARCH FUNCTION*/

void search(float nSCH, float sigmaSCH, int CxSCH, int CySCH, float rminSCH, float rmaxSCH, int respSCH, int resiSCH){
int arrszSCH=rmimxPD(rminSCH, rmaxSCH);
rmimxbdPD(rminSCH, rmaxSCH, arrszSCH);
memset(maxradSCH, 0, sizeof(maxradSCH));
memset(maxblurSCH, 0, sizeof(maxblurSCH));
/*Making a 10x10 neighbourhood around the ordinate input CxSH and CySH*/
int iteri, iterj;
for (iteri=(CxSCH-5);iteri<(CxSCH+6);iteri++) {
for (iterj=(CySCH-5);iterj<(CySCH+6);iterj++) {
		int idx=partiald(iteri, iterj, rminSCH, rmaxSCH, sigmaSCH, nSCH, respSCH, resiSCH);
		maxradSCH[iteri][iterj]=*(R+idx);
		maxblurSCH[iteri][iterj]=*(resmat+idx);
		/*printf("\niteri:%d  iterj:%d  idx:%d\n",iteri,iterj,idx);
		printf("\nR_PD:%f   Blur_PD:%f  iteri:%d  iterj:%d\n",maxradSCH[iteri][iterj],maxblurSCH[iteri][iterj],iteri,iterj);*/
}
}
*(FINE_SCH+2)=max2d(maxblurSCH);
}

/*THRESH FUNCTION*/

void thresh(float thcol, float rminTHS, float rmaxTHS, float sigmaTHS,
			float nTHS){
int iteri, iterj, size_thcol=0;
for (iteri=0;iteri<ROWS;iteri++) {
for (iterj=0;iterj<COLS;iterj++) {
if ((array[iteri][iterj])<thcol){size_thcol=size_thcol+1;}
}
}
X_th=(int*)malloc(size_thcol*sizeof(int));
Y_th=(int*)malloc(size_thcol*sizeof(int));
int *tempkeep;
tempkeep=(int*)malloc(size_thcol*sizeof(int));
int seliter=0;
for (iteri=0;iteri<ROWS;iteri++) {
for (iterj=0;iterj<COLS;iterj++) {
if ((array[iteri][iterj])<thcol){*(X_th+seliter)=iteri;*(Y_th+seliter)=iterj;seliter=seliter+1;
}
}
}
int pix_sel; int net_pxlcnt=0;
for (pix_sel=0;pix_sel<size_thcol;pix_sel++) {
	if (((*(X_th+pix_sel))>(1.3*rminTHS))&&((*(Y_th+pix_sel))>(1.3*rminTHS))&&
	((*(X_th+pix_sel))<=(1.3*(ROWS-rminTHS)))&&((*(Y_th+pix_sel))<=(1.3*(COLS-rminTHS))))
	{
		int i,j;
		float mymin = array[*(X_th+pix_sel)][*(Y_th+pix_sel)];
			for(i = (*(X_th+pix_sel))-1; i < (*(X_th+pix_sel))+2; i++)	{
				for(j = (*(Y_th+pix_sel))-1; j < (*(Y_th+pix_sel))+2; j++) {
					if(array[i][j] < mymin)
					{
						mymin = array[i][j];
					}
				}
			}
		if ((array[*(X_th+pix_sel)][*(Y_th+pix_sel)])==mymin) {
		*(tempkeep+net_pxlcnt)=pix_sel; net_pxlcnt=net_pxlcnt+1;
		}
	}
}
X_fn=(int*)malloc(net_pxlcnt*sizeof(int));
Y_fn=(int*)malloc(net_pxlcnt*sizeof(int));
int iter_fn;
for (iter_fn=0;iter_fn<net_pxlcnt;iter_fn++) {
*(X_fn+iter_fn)=*(X_th+(*(tempkeep+iter_fn)));
*(Y_fn+iter_fn)=*(Y_th+(*(tempkeep+iter_fn)));
}
memset(maxradSCH, 0, sizeof(maxradSCH));
memset(maxblurSCH, 0, sizeof(maxblurSCH));
free(tempkeep);free(X_th);free(Y_th);
for (iter_fn=0;iter_fn<net_pxlcnt;iter_fn++) {
int pdret=partiald((*(X_fn+iter_fn)), (*(Y_fn+iter_fn)), rminTHS,
					rmaxTHS, 300, nTHS, 0, 1);
/*I am setting the value of sigma to 300*/
maxblurSCH[*(X_fn+iter_fn)][*(Y_fn+iter_fn)]=*(resmat+pdret);
maxradSCH[*(X_fn+iter_fn)][*(Y_fn+iter_fn)]=*(R+pdret);
}

/* Code for printing (DISABLED!! */
       /* int i3,j3;
        for ( i3 = 0; i3 < ROWS; i3++ )
      {
         puts("\n");
		 printf("\narray[%lu]: ", (long unsigned)i3);
		 puts("\n");
         for ( j3 = 0; j3 < COLS; j3++ )
         {
            printf("%f ", maxblurSCH[i3][j3]);
         }
		}*/
float max_ret=max2d(maxblurSCH);
int X_1=*(FINE_SCH+0);int Y_1=*(FINE_SCH+1);
search(nTHS, sigmaTHS, X_1, Y_1, rminTHS, rmaxTHS, 0, 1);
/* 0 = pupil "yes" :: 1 = iris "no" and vice versa */
int X_1pup=(int)*(FINE_SCH+0);int Y_1pup=(int)*(FINE_SCH+1);
float RAD_pup=*(FINE_SCH+2);
search(nTHS, sigmaTHS, X_1pup, Y_1pup, 1.9*RAD_pup, 5*RAD_pup, 1, 0);
int X_1iris=(int)*(FINE_SCH+0);int Y_1iris=(int)*(FINE_SCH+1);
float RAD_iris=*(FINE_SCH+2);
printf("\nThe values obtained are:\nXpup=%d, Ypup=%d, Rpup=%f\n",X_1pup, Y_1pup, RAD_pup);
printf("\nThe values obtained are:\nXiris=%d, Yiris=%d, Riris=%f\n",X_1iris, Y_1iris, RAD_iris);
}

/*START OF THE MAIN FUNCTION*/

int main(void) {
arraybuild();
FINE_SCH=(float*)malloc(3*sizeof(float));
/*VALUE INITIALISATION PART*/
float n=600; /*float r;*/
/*int C[2]={0}; char part[8]={0};    /*For storing the input to search for*/
/*printf("Enter the value of n\n"); scanf("%f", &n);*/
/*printf("Enter the value of r\n"); scanf("%f", &r);*/
/*printf("Enter the part to scan for\n"); scanf("%s", &part);*/
/*printf("The part you entered is    "); puts(part);*/
/*char pupil[]="pupil"; char iris[]="iris";  /*For comparison*/
/*int respupil=strcmp(part,pupil);  /*Immediately storing the result of comparisons*/
/*int resiris=strcmp(part,iris);    /*Immediately storing the result of comparisons*/
/*printf("Enter the x coordinate\n"); scanf("%d", &C[0]);*/
/*printf("Enter the y coordinate\n"); scanf("%d", &C[1]);*/
/*printf("Enter the threshold value: for selecting black pixels(0.5 or 127)\n");*/
float pixintcol=127;
thresh(pixintcol, 16, 30, 0.5, n);
return 0;
}
