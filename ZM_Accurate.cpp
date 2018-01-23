/********************************************************************************************
*             Accurate calculation of Zernike Moments(ZMs).   								*
*___________________________________________________________________________________________*
*	Developed By: Dr. Chandan Singh and Jaspreet Singh.										*
*             Dr. Chandan Singh																*
*		      Professor(Re-employed),														*
*		      Department of Computer Science,												*
*		      Punjabi University, Patiala, India - 147 002									*
*             Mobile : +91-9872043209, email : chandan.csp@gmail.com                    	*
*             Jaspreet Singh																*
*			  Research Scholar (Ph.D),                                                      *
*             Department of Computer Science,                                               *
*			  Punjabi University, Patiala, India - 147 002                                  *
*			  Mobile : +91-9803335164, email : jaspreetmaancs@gmail.com                 *
*     	Date : September 21,2017.															*
*	Reference : 1) C. Singh, E. Walia, R. Upneja, Accurate calculation of Zernike moments,  *
*				   Inf. Sci. (Ny). 233 (2013) 255–275.                                      *
*				2) TIFF, Revision 6.0, (June 1992). Adobe Developers Association.			*
*	               www.adobe.com/Support/TechNotes.html										*
*                  and ftp://ftp.adobe.com/pub/adobe/DeveloperSupport/TechNotes/PDFfiles	*
*               3) Digital Image Processing by Rafael C. Ganzalez and						*
*                  Richard E. Woods, Third Edition, Pearson Education, India, (for image	*
*                  processing operations)													*
*___________________________________________________________________________________________*
*               The program reads the TIFF file,interprets,processes and writes				*
*               the processed image to an output file.										*
*___________________________________________________________________________________________*
*      Notes: 1) Program has been implemented in Microsoft's Visual C++6.0 under			*
*                Windows environment and can be execetued under "Build" menu				*
*                item of Visual C++6.0 or through "command line". To provide				*
*                sufficient memory for its compilation and execution, under Integrated      *
*                Development Environment(IDE) choose the options:                           * 
*				 Projects->Settings->Link->Category->output->Reserve(insert the value,      *
*                say, 1000000000 for one GB of Stack). Without this setting                 *
*                you may fail to execute the program. Also, it may not run					*
*                purely under DOS environment because the memory requirement is				*
*                high.																		*
*             2) This is not a product.														*
*___________________________________________________________________________________________*
*********************************************************************************************
*/
#include <iostream>
#include <stdio.h>
#include <conio.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#define MAXSIZE 512 // Maximum size of the image
#define NMAX 60  // Maximum number of ZMs
const double EPS=1.0e-06;  // Tolerence, used in the computation of PSNR
using namespace std;

void  Read_Tiff_File(void);
void copy_array_to_file(unsigned long ImageLength, unsigned long ImageWidth, 
				   unsigned char fxyout[MAXSIZE][MAXSIZE], FILE *outfile_ptr);
void copy_to_array(unsigned long ImageLength, unsigned long ImageWidth, 
				   unsigned char *IB, unsigned char fxy[MAXSIZE][MAXSIZE]);
void display_error( char *message);
void Zernike_moment_numerical_integration(int N, int nmax, unsigned char fxy[MAXSIZE][MAXSIZE],
					double Anm[NMAX+1][2*NMAX+1], double AnmReal[NMAX+1][2*NMAX+1], 
					double AnmImg[NMAX+1][2*NMAX+1], FILE *momptr, unsigned char in_filename[]);
void Zernike_moment_symmetry8(int N, int nmax, unsigned char fxy[MAXSIZE][MAXSIZE],
					double Anm[NMAX+1][2*NMAX+1], double AnmReal[NMAX+1][2*NMAX+1], 
					double AnmImg[NMAX+1][2*NMAX+1], FILE *momptr );
void Zernike_moment_accurate(int N, int pmax, unsigned char fxy[MAXSIZE][MAXSIZE],
					double Anm[NMAX+1][2*NMAX+1], double AnmReal[NMAX+1][2*NMAX+1], 
					double AnmImg[NMAX+1][2*NMAX+1], FILE *momptr, unsigned char in_filename[]);
void Reconstruction_ZM(int N, int nmax, double Anm[NMAX+1][2*NMAX+1], double AnmReal[NMAX+1][2*NMAX+1], 
					double AnmImg[NMAX+1][2*NMAX+1], unsigned char fxyoriginal[MAXSIZE][MAXSIZE],
					unsigned char fxycouter[MAXSIZE][MAXSIZE], FILE *momptr, unsigned char in_filename[14]);
void MSRE(int N, unsigned char fxyoriginal[MAXSIZE][MAXSIZE], double fxyR[MAXSIZE][MAXSIZE], FILE *momptr);

unsigned long ImageWidth, ImageLength, StripOffsets, RowsPerStrip,
StripsPerImage, ImageSize, StripByteCounts;
double Anm[NMAX+1][2*NMAX+1]={0}, AnmReal[NMAX+1][2*NMAX+1]={0}, AnmImg[NMAX+1][2*NMAX+1]={0}; 
double xc, yc, Rad, xbarf, ybarf, Df, Dfsqr, cpu_time_used;
const double pi=4.0*atan(1.0), pi2=2*pi;
const int LMAX=MAXSIZE*MAXSIZE*(NMAX+1)*(NMAX+1)/2;
const int MAXINT=20;  // maximum no. of Gaussian integtration points
unsigned char **new_image;
short int BitWhite;
short int ResolutionUnit, BitsPerSample;
unsigned long XResolution, YResolution, NoOfStrips;
void *ImageBuffer;
clock_t start, stop;

int main()
{
	Read_Tiff_File();
	return(0);
}

void  Read_Tiff_File()
{
	long int n, nmax;
	unsigned int i, j, reply;
	short unsigned int byte_order, TIFF_identifier;
	unsigned short int no_of_entries, tag, field_type, temp;
	unsigned long offset_address, count, value, so, sbc, StripAddress,
		StripSize;
	unsigned long file_position, N, N2;
	unsigned char *IB, item;
	unsigned char *ib, fxy[MAXSIZE][MAXSIZE], fxyout[MAXSIZE][MAXSIZE];
	// Reads an image file in the variable "in_filename" and writes the output to 
	// the output image file "out.tif".
	unsigned char in_filename[25] = "lena128.tif", out_filename[25]="out.tif"; 
	unsigned char fxyoriginal[MAXSIZE][MAXSIZE];
	FILE *infile_ptr, *outfile_ptr, *momptr;
	
	// Open the input file for processing
	// printf("Input TIFF File ? : ");
	// scanf( "%s", in_filename);
	if((infile_ptr = fopen((const char *)in_filename, "rb+")) ==NULL)
	{
		printf("can't open file - %s - either the filename is wrong or the"
			" file does not exist\n", in_filename);
		getch(); exit(0);
	}
	// Open a new file for storing the processed image
	// printf("Output TIFF File ? : ");
	// scanf( "%s", out_filename);
	if((outfile_ptr = fopen((const char *)out_filename, "wb")) ==NULL)
	{
		printf("Error in opening file - %s - the filename may be wrong"	, out_filename);
		getch(); exit(0);
	}
	// The file "moment.out" is used for saving the ZMs.
	if((momptr = fopen("moment.out", "w")) ==NULL) 
	{
		printf("Error in opening file - moment.out");
		getch(); exit(0);
	}
	// Copy the existing file to the new file  
	while( fread((void *)&item, (size_t)1, (size_t)1, infile_ptr) !=0)
	{
		fwrite((void *)&item, (size_t)1, (size_t)1, outfile_ptr);
	}
	// Bring the file pointer to the beginning of the input file
	rewind(infile_ptr);
	// Close the output file and reopen it in read & write modes
	fclose(outfile_ptr);
	if((outfile_ptr = fopen((const char *)out_filename, "rb+")) ==NULL)
	{
		printf("Error in opening file - %s - the filename may be wrong"
			, out_filename);
		getch(); exit(0);
	}
	// printf("Enter 1 for Zernike Moments Using ZOA 8-Way Symmetry\nEnter 2 for Accurate"
	//    " Computation of ZM Using GM\n"
	//	"Enter 3 for Accurate Computation of ZM Using Numerical Integration:");
	// scanf("%d", &reply); 
	// reply=1 for computation of Zernike Moments Using zeroth order approximation with 
	// 8-Way Symmetry.
	// reply=2 for Accurate Computation of ZM Using GM.
	// reply=3 for Accurate Computation of ZM Using Numerical Integration.
	reply=1;
	if(reply<1 || reply>3) display_error("Choice no. is wrong.");
	
	fread((void *)&byte_order, (size_t)2, (size_t)1, infile_ptr);
	fread((void *)&TIFF_identifier, (size_t)2, (size_t)1, infile_ptr);
	if(byte_order != 0x4949)  display_error("Byte Order not Intel Processor"
		" Compatible");
	if(TIFF_identifier != 0x002A)  display_error("File is not a valid TIFF file");
	fread((void *)&offset_address, (size_t)4, (size_t)1, infile_ptr);
	fseek(infile_ptr, offset_address, SEEK_SET);
	
	fread((void *)&no_of_entries, (size_t)2, (size_t)1, infile_ptr);
	
	for(i=0; i<no_of_entries; i++)
	{
		fread((void *)&tag, (size_t)2, (size_t)1, infile_ptr);
		fread((void *)&field_type, (size_t)2, (size_t)1, infile_ptr);
		fread((void *)&count, (size_t)4, (size_t)1, infile_ptr);
		fread((void *)&value, (size_t)4, (size_t)1, infile_ptr);
		switch(tag)
		{
		case 256 : 
			ImageWidth = value;
			break;
		case 257 : 
			ImageLength = value;
			break;
		case 258 :
			BitsPerSample = (short)value;
			break;
		case 259 : 
			if(value !=1) display_error("Image is compressed. "
			"Can't process");
			break;
			case 262 : 
			if(value<0 || value >3)
			display_error("The Image is not binary/8-bit gray"
				"/8-bit coloured/24-bit colour");
			BitWhite =0;
			if(value ==1) BitWhite =1;
			break;
					
		case 273 : 
			NoOfStrips = count;
			StripOffsets= value;
			break;
		case 278 : 
			RowsPerStrip = value;
			StripsPerImage =(ImageLength + RowsPerStrip -1)/RowsPerStrip;
			break;
		case 279 : 
			StripByteCounts = value;
			break;
		case 282 : 
			XResolution = (short)value;
			break;
		case 283 : 
		    YResolution = (short)value;
			break;
		case 296 : 
			ResolutionUnit= (short)value;
			break;
		default :   break;
		}
	}
	
	if(BitsPerSample !=8) display_error("Image is not a gray scale image");
	ImageSize = ImageWidth*ImageLength;
	n=ImageSize;
	// printf("ImageSize=%lu\n", ImageSize);
	ImageBuffer =  malloc(ImageSize);
	if(ImageBuffer == NULL)
		display_error("Unable to allocate memory for storing the image"
		" in memory");
	
	// Zeroize the contents of the Image Buffer
	IB = (unsigned char *)ImageBuffer;
	for(i=0; i<ImageLength; i++)
		for(j=0; j<ImageWidth; j++)
		{
			(*IB) =0;
			IB++;
		}
		
		so = StripOffsets;  sbc = StripByteCounts;  ib = (unsigned char *)ImageBuffer;
		if(NoOfStrips ==1)
		{
			fseek(infile_ptr, so, SEEK_SET);
			file_position = (unsigned long)ftell(infile_ptr);
			temp =fread(ImageBuffer, (size_t)1, (size_t)ImageSize, infile_ptr);
		}
		else
		{
			for(i=0; i<NoOfStrips; i++)
			{
				fseek(infile_ptr, so, SEEK_SET);
				fread(&StripAddress, 4, 1, infile_ptr);
				so = so + 4;
				fseek(infile_ptr, sbc, SEEK_SET);
				sbc = sbc +4;
				fread(&StripSize, 4, 1, infile_ptr);
				
				fseek(infile_ptr, StripAddress, SEEK_SET);
				if(i==0)
				{
					file_position = (unsigned long) ftell(infile_ptr);
				}
				fread(ib, (size_t)StripSize, (size_t)1, infile_ptr);
				(unsigned char *)ib = (unsigned char *)ib + StripSize;
			}
		}
		IB = (unsigned char *)ImageBuffer;
        copy_to_array(ImageLength, ImageWidth, IB, fxy);
		// printf("Enter the order of moments:");
		// scanf("%d",&nmax);
		if( ImageLength != ImageWidth) display_error("Image is not a square image");
		N = ImageLength;
		if( N>MAXSIZE ) display_error("Size of image is more than allowed");
		if( N%2==1 ) printf("Image size is odd\n");
		nmax=40; // order of moments
		if( nmax > 60 )
		{
			display_error("Entered order of moments exceeds the maximum limit");
		}
		if( reply==1 )
		{
			Zernike_moment_symmetry8(N, nmax, fxy, Anm, AnmReal, AnmImg, momptr);
			Reconstruction_ZM(N, nmax, Anm, AnmReal, AnmImg, fxy, fxyout,  momptr, 
				in_filename);
		}
		if( reply==2 )
		{
			Zernike_moment_accurate(N, nmax, fxy, Anm, AnmReal, AnmImg, momptr, in_filename );
			Reconstruction_ZM(N, nmax, Anm, AnmReal, AnmImg, fxy, fxyout,  momptr, in_filename);
		}
		if( reply==3 )
		{
			Zernike_moment_numerical_integration(N, nmax, fxy, Anm, AnmReal, AnmImg, momptr, in_filename);
			Reconstruction_ZM(N, nmax, Anm, AnmReal, AnmImg, fxy, fxyout, momptr, in_filename);
		}		
		fseek(outfile_ptr, file_position, SEEK_SET);
        copy_array_to_file(ImageLength, ImageWidth, fxyout, outfile_ptr);
		fclose(infile_ptr);  fclose(outfile_ptr);
        printf("Press a key to exit the program\n");		
		getch();  exit(0);
		// Operation on the input image file is over
}

/*
/*************************************************************************
*                                                                        *
*   Calculation of Zernike moments using 8-way Symmetry property.        *
*   q-recursion is used to compute radial polynomials.                   *
**************************************************************************
*/
void Zernike_moment_symmetry8(int N, int nmax, unsigned char fxy[MAXSIZE][MAXSIZE],
					double Anm[NMAX+1][2*NMAX+1], double AnmReal[NMAX+1][2*NMAX+1], 
					double AnmImg[NMAX+1][2*NMAX+1], FILE *momptr)
{
	int n, m, x, y, icount, xbar, ybar, x1, x2, x3, x4, y1, y2, y3,
		y4, RAD, f0, f1, f2, f3, f4, f5, f6, f7, option=1;
	double mult, xf, yf, rsqr, rsqrt, pi=4.0*atan(1.0), area, Rn, Rnm, 
		Rnm2, Rnmp2, Rnmp4, Rnnm2, H1, H2, H3, np1byarea, a, b, C, S, 
		temp, cost, sint, N2, rpower[NMAX+1], H1mn[NMAX+1][NMAX+1],  
		H2mn[NMAX+1][NMAX+1], H3mn[NMAX+1][NMAX+1], sinmt[NMAX+1], cosmt[NMAX+1];
	
	RAD=N/2;  	
	N2=N*sqrt(2.0);  // outer circle, for inner circle use N2=N
	fprintf(momptr, "\nOuter circle is used\n");

	Dfsqr=(N2/2.0)*(N2/2.0);
   	area=pi*Dfsqr; 
	//printf("Calculation of moments using zeroth order approximation and 8-way symmetry\n");
	fprintf(momptr, "\nCalculation of moments using zeroth order approximation and 8-way symmetry");
	fprintf(momptr, "\nOuter circle is used\n");
	start=clock(); // start clock time
	// Initialize the Zernike moments.
	for(n=0; n<=nmax; n++)
		for(m=0; m<=nmax; m++)
		{
			AnmReal[n][m]=0.0; AnmImg[n][m]=0.0;
			Anm[n][m]=0.0;
		}
    xbar=(int)RAD;  ybar=(int)RAD;
	// Compute H1, H2 and H3 and save. Values for m=n and m=n-2 are skipped
	for(n=0; n<=nmax; n++)
	{
		for(m=n-4; m>=0; m=m-2)
		{
			H3mn[n][m]=-(double)(4*(m+2)*(m+1))/((n+m+2)*(n-m));
			H2mn[n][m]=H3mn[n][m]*(n+m+4)*(n-m-2)/(4*(m+3))+(m+2);
			H1mn[n][m]=(m+4)*(m+3)/(double)2-H2mn[n][m]*(m+4)+
				H3mn[n][m]*(n+m+6)*(n-m-4)/8.0;
		}  // for m-loop
	}  // for n-loop
	// If the centre of the circle is (h,k), then 8-way symmetry gives the points :
	// 1-(h+x, k+y), 2-(h+y, k+x), 3-(h-y, k+x), 4-(h-x, k+y), 5-(h-x, k-y), 
	// 6-(h-y, k-x), 7-(h+y, k-x), 8-(h+x, k-y), 
	for(x=0; x<RAD; x++)
	{
		xf=(2.0*x+1)/N2;
		x1=xbar+x;  x2=xbar-x-1;  y3=ybar+x;  
		y4=ybar-x-1;
		for(y=0; y<=x; y++)
		{
			yf=(2.0*y+1)/N2;
			rsqr=xf*xf+yf*yf;
			if(rsqr>1.0) continue;
			rsqrt=sqrt(rsqr);
			y1=ybar+y;  y2=ybar-y-1;  x3=xbar+y;  
			x4=xbar-y-1;
			// Skip if any one of the indices is outside range
			if(x1>=N || x2<0 || x3>=N || x4<0 || 
				y1>=N || y2<0 || y3>=N || y4<0) continue;  
			f0=fxy[x1][y1];  f1=fxy[x3][y3];  f2=fxy[x4][y3]; 
			f3=fxy[x2][y1];  f4=fxy[x2][y2];  f5=fxy[x4][y4];
			f6=fxy[x3][y4];  f7=fxy[x1][y2];
			// Computue r**n
			mult=rsqrt;
			rpower[0]=1.0;
			for(n=1; n<=nmax; n++)
			{
				rpower[n]=mult;
				mult=mult*rsqrt;
			}  // for n-loop
			// Compute cos(m*theta) and sin(m*theta) for m=0
			cosmt[0]=1.0;  	sinmt[0]=0.0;
            // for m=1
			if(rsqrt==0.0) { a=1.0;  b=0.0; }  // theta is 0.0 for xf=0.0 & yf=0.0
			else { a=xf/rsqrt;  b=yf/rsqrt; }
			cosmt[1]=a;  
			sinmt[1]=b;
			C=a;  S=b;
			// for m>=2
			for(m=2; m<=nmax; m++)
			{
				temp=a*C-b*S;
				S=a*S+b*C;
				cosmt[m]=temp;  
				sinmt[m]=S;
				C=temp;  
			}  // m-loop
			for(n=0; n<=nmax; n++)
			{
				Rn=rpower[n];
				if(n>1) Rnm2=rpower[n-2];
				for(m=n; m>=0; m=m-2)
				{
					if(m==n)
					{
						Rnm=Rn;
						Rnmp4=Rn;  // save Rn(m+4)
					}
					else if(m==n-2)
					{
						Rnnm2=n*Rn-(n-1)*Rnm2;
						Rnm=Rnnm2;
						Rnmp2=Rnnm2;  // save Rn(m+2)
					}
					else
					{
						H3=H3mn[n][m];
						H2=H2mn[n][m];
						H1=H1mn[n][m];
						if(rsqr==0.0)
						{
							if(m==0)
							{
								if(n%4==0) Rnm=1.0;
								else Rnm=-1.0;
							}
							else Rnm=0.0;
						}  // if-rsqr
						else Rnm=H1*Rnmp4+(H2+H3/rsqr)*Rnmp2;
						Rnmp4=Rnmp2;  // Rnmp2 now becomes Rnmp4 for next m
						Rnmp2=Rnm;  // Rnm becomes Rnmp2 for next m
					}  // if block
					cost=cosmt[m];  sint=sinmt[m];
					if(N%2==1 && x==0)  // for(0,0) which for odd size
					{
						AnmReal[n][m]=AnmReal[n][m]+Rnm*f0;
					}  // x==0
					else if(x==y)  // for 45 deg
					{
						switch(m%4)
						{
						case 0 :
							AnmReal[n][m]=AnmReal[n][m]+(f0+f2+f4+f6)*cost*Rnm;
							AnmImg[n][m]=AnmImg[n][m]+(f0+f2+f4+f6)*sint*Rnm;
							break;
						case 1 :
							AnmReal[n][m]=AnmReal[n][m]+((f0-f4)*cost+(-f2+f6)*sint)*Rnm;
							AnmImg[n][m]=AnmImg[n][m]+((f2-f6)*cost+(f0-f4)*sint)*Rnm;
							break;
						case 2 :
							AnmReal[n][m]=AnmReal[n][m]+(f0-f2+f4-f6)*cost*Rnm;
							AnmImg[n][m]=AnmImg[n][m]+(f0-f2+f4-f6)*sint*Rnm;
							break;
						case 3 :
							AnmReal[n][m]=AnmReal[n][m]+((f0-f4)*cost+(f2-f6)*sint)*Rnm;
							AnmImg[n][m]=AnmImg[n][m]+((-f2+f6)*cost+(f0-f4)*sint)*Rnm;
							break;
						}  // switch
					}
					else   // for general case
					{
						switch(m%4)
						{
						case 0 :
							AnmReal[n][m]=AnmReal[n][m]+(f0+f1+f2+f3+f4+f5+f6+f7)*Rnm*cost;
							AnmImg[n][m]=AnmImg[n][m]+(f0-f1+f2-f3+f4-f5+f6-f7)*Rnm*sint;
							break;
						case 1 :
							AnmReal[n][m]=AnmReal[n][m]+((f0-f3-f4+f7)*cost+
									            (f1-f2-f5+f6)*sint)*Rnm;
							AnmImg[n][m]=AnmImg[n][m]+((f1+f2-f5-f6)*cost+
								(f0+f3-f4-f7)*sint)*Rnm;
							break;
						case 2 :
							AnmReal[n][m]=AnmReal[n][m]+(f0-f1-f2+f3+f4-f5-f6+f7)*Rnm*cost;
							AnmImg[n][m]=AnmImg[n][m]+(f0+f1-f2-f3+f4+f5-f6-f7)*Rnm*sint;
							break;
						case 3 :
							AnmReal[n][m]=AnmReal[n][m]+((f0-f3-f4+f7)*cost+
									            (-f1+f2+f5-f6)*sint)*Rnm;
							AnmImg[n][m]=AnmImg[n][m]+((-f1-f2+f5+f6)*cost+
								(f0+f3-f4-f7)*sint)*Rnm;
							break;
						}  // switch
					}  // else of general case
				}  // for m-loop
			}  // for n-loop
		}  // for y-loop
	}  // for x-loop
	// Compute ZMs.
	stop = clock(); // stop clock time
	fprintf(momptr, "\n    n    m    AnmReal     AnmImg       Anm        \n");
	for(n=0; n<=nmax; n++)
	{
		np1byarea=(n+1)/area;
		icount=0;
		for(m=0; m<=n; m++)
		{
			// n-|m| must be even for ZMs.
			if((n-abs(m))%2 !=0) continue; 
			AnmReal[n][m]=np1byarea*AnmReal[n][m];
			AnmImg[n][m]=-np1byarea*AnmImg[n][m];
			Anm[n][m] = sqrt(AnmReal[n][m]*AnmReal[n][m]+AnmImg[n][m]*AnmImg[n][m]);
			if(m>=0)
			{
				if(icount==0)
				{
					fprintf(momptr, "%5d%5d%12.6f%12.6f%12.6f\n", 
						n, m, AnmReal[n][m], AnmImg[n][m], Anm[n][m]); 
					icount++;
				}
				else
					fprintf(momptr, "%10d%12.6f%12.6f%12.6f\n", 
					m, AnmReal[n][m], AnmImg[n][m], Anm[n][m]); 
			}  
		} // next m
	}  // next n
	cpu_time_used = ((double) (stop - start)) / CLOCKS_PER_SEC;
	fprintf(momptr, "Time taken by Zernike_moment_symmetry8() for the computation"
		" of moments is : %f sec\n", cpu_time_used); 
	return;
}

/*
/*************************************************************************
*                                                                        *
*   Numerical double integration using Gaussian quadrature.              *
*                                                                        *
**************************************************************************
*/
double integration(double a, double b, double c, double d, double n, double m)
{
	int i, j, NW, NW2;
	double p, q, r, t, xi, xsqr, yj, integral, fvalue, xvalue, wi, 
		w[100], x[100];
	NW=5;  NW2=NW/2;
	x[0]=-0.90617985;  x[1]=-0.53846931;  x[2]=0.0;
	w[0]=0.23692688;   w[1]=0.47862867;   w[2]=0.56888889;
	
	for(i=1; i<=NW2; i++)
	{
		x[NW2+i]=-x[NW2-i];  w[NW2+i]=w[NW2-i];
	}
	p=(a+b)/2.0;  q=(b-a)/2.0;  r=(c+d)/2.0;  t=(d-c)/2.0;
	integral=0.0;
	for(i=0; i<NW; i++)
	{
		xi=p+q*x[i];  xsqr=xi*xi;  xvalue=pow(xi, n);  wi=w[i];  
		for(j=0; j<NW; j++)
		{
			yj=r+t*x[j];
			fvalue=xvalue*pow(yj, m);
            integral=integral+wi*w[j]*fvalue;
		}  // for j-loop
	}  // for i-loop
	integral=q*t*integral;
	return integral;
}

/*
/*************************************************************************
*                                                                        *
*   Exact calculation of exact Zernike Moments using Geomentric moments.  *
*   Geometric moments are also computed using numerical integration just *
*   for the sake of verifying the accuracy of numerical integration. The *
*   gray values are assumed to be constant over a pixel grid.            *
*                                                                        *
**************************************************************************
*/
void Zernike_moment_accurate(int N, int pmax, unsigned char fxy[MAXSIZE][MAXSIZE],
					double Anm[NMAX+1][2*NMAX+1], double AnmReal[NMAX+1][2*NMAX+1], 
					double AnmImg[NMAX+1][2*NMAX+1], FILE *momptr, unsigned char in_filename[25] )
{
	int i, j, k, m, n, p, q, s, icount, pmkby2, ppkby2, kmqby2, kpqby2, 
		m2n, option=1, option_integration=1; 
	double D, start, xi, del, delby2, xipdelby2, ximdelby2, 
		yj, yjpdelby2, yjmdelby2, sum, area, pp1byarea, 
		sumR, sumI, bpqk, term, term1, term2, hpi, 
		factorial[2*NMAX+1], xu[MAXSIZE], xl[MAXSIZE],
		B[NMAX+1][NMAX+1][NMAX+1], H[NMAX+1][MAXSIZE], M[NMAX+1][NMAX+1], 
		sR[NMAX+1][NMAX+1], sI[NMAX+1][NMAX+1];
	double pi=4.0*atan(1.0);
	printf("Exact Calculation of Zernike Moments using Geomentric moments\n");
	fprintf(momptr, "Exact Calculation of Zernike Moments using Geomentric moments\n"
		"Note:Outer circle is used\n");
	D=N*sqrt(2.0);  // D=N for inner circle
    area=pi;
	// printf("Enter 1 for Computation of Geomentric moments without symmetry"
	//	"\nEnter 2 for Computation of Geometric Moments using numerical integration:");
	// scanf("%d",&option_integration);
	// option_integration = 1 for Computation of Geomentric moments without symmetry
	// option_integration = 2 for Computation of Geometric Moments using numerical integration
	option_integration = 2;
	start=clock(); // start clock time
	// Initialize the Zernike moments.
	for(p=0; p<=pmax; p++)
		for(q=0; q<=pmax; q++)
		{
			AnmReal[p][q]=0.0; AnmImg[p][q]=0.0;
			Anm[p][q]=0.0;
		}
    // Step-I: Compute all factorial values once
	factorial[0]=1.0;for(i=1; i<=2*pmax+1; i++) factorial[i]=i*factorial[i-1];
	// Step-II: Compute B[p][q][k]
	for(p=0; p<=pmax; p++)
	{
		for(q=0; q<=pmax; q++)
		{
			for(k=q; k<=p; k++)
			{
				pmkby2=(p-k)/2;   ppkby2=(p+k)/2;
				kmqby2=(k-q)/2;   kpqby2=(k+q)/2;
				term=factorial[ppkby2]/
					(factorial[pmkby2]*factorial[kpqby2]*factorial[kmqby2]);
				if(pmkby2%2==0) B[p][q][k]=term;
				else B[p][q][k]=-term;
			} // k-loop
		}  // q-loop
	}  // p-loop

	// Step-III: Compute H[p][x] or H[q][y] : both are same for square image.
	// Outer circle is considered.
	del=2.0/D;  delby2=del/2.0;  start=-(1+N)/D; 
	for(i=0; i<N; i++)
	{
		xi=(2*i+1-N)/D;
		xipdelby2=xi+delby2;  ximdelby2=xi-delby2;
		xu[0]=xipdelby2;  xl[0]=ximdelby2;  // for p=0, p+1=1
		for(p=1; p<=pmax+1; p++)
		{
			xu[p]=xu[p-1]*xipdelby2;  xl[p]=xl[p-1]*ximdelby2;
		}
		for(p=0; p<=pmax; p++)
		{
			H[p][i]=(xu[p]-xl[p])/(p+1);
		}  // i-loop
	}  // p-loop
		
	// Step-IV: Computation of Geometric Moments, M[p][q]
	if(option_integration==1)
	{
		// Computation of Geomentric moments, M[p][q], without symmetry
		for(p=0; p<=pmax; p++)
		{
			for(q=0; q<=pmax; q++)
			{
				sum=0.0;
				for(i=0; i<N; i++)
				{
					xi=(2*i+1-N)/D;  hpi=H[p][i];
					for(j=0; j<N; j++)
					{
						yj=(2*j+1-N)/D;
						sum=sum+fxy[i][j]*hpi*H[q][j];
					}  // i-loop
				}  // j-loop
				// printf(" without symmetry: p=%d, q=%d, sum=%f\n", p, q, sum);  getch();
				M[p][q]=sum;
				// printf(" without symmetry: p=%d, q=%d, M[p][q]=%f\n", p, q, M[p][q]);  
				// getch();
			}  // q-loop
		}  // p-loop
	}
	else if(option_integration==2)
	{
		// Step-IV: Computation of Geometric Moments, M[p][q], 
		// using numerical integration
		printf("Using numerical integration : Wait as it takes time\n");  
		for(p=0; p<=pmax; p++)
		{
			for(q=0; q<=pmax; q++)
			{
				sum=0.0;
				yj=start;  
				for(j=0; j<N; j++)
				{
					yj=yj+del;
					yjpdelby2=yj+delby2;  yjmdelby2=yj-delby2;
					xi=start;
					for(i=0; i<N; i++)
					{
						xi=xi+del;
						xipdelby2=xi+delby2;  ximdelby2=xi-delby2;
						sum=sum+fxy[j][i]*integration(ximdelby2, xipdelby2, 
							yjmdelby2, yjpdelby2, (double)p, (double)q);
					}  // i-loop
				}  // j-loop
				M[p][q]=sum;
			}  // q-loop
		}  // p-loop
	}
    // Step-V: Computation of intermediate terms sR[q][k] and sI[q][k]
	for(q=0; q<=pmax; q++)
	{
		for(k=q; k<=pmax; k++)
		{
				s=(k-q)/2;
				sumR=0.0;  sumI=0.0;
				for(m=0; m<=s; m++)
				{
					term1=factorial[s]/(factorial[m]*factorial[s-m]);
					for(n=0; n<=q; n++)
					{
						m2n=2*m+n;
						term2=factorial[q]/(factorial[n]*factorial[q-n]);
						term=term1*term2*M[k-m2n][m2n];
						if(q<=0)
						{
							switch(n%4)
							{
							case 0 : sumR=sumR+term;  break;
							case 1 : sumI=sumI+term;  break;
							case 2 : sumR=sumR-term;  break;
							case 3 : sumI=sumI-term;  break;
							} // switch
						} // if
						else
						{
							switch(n%4)
							{
							case 0 : sumR=sumR+term;  break;
							case 1 : sumI=sumI-term;  break;
							case 2 : sumR=sumR-term;  break;
							case 3 : sumI=sumI+term;  break;
							} // switch
						} // else
					} // n-loop
				} // m-loop
				sR[q][k]=sumR;  sI[q][k]=sumI;
		}  // k-loop
	}  // q-loop
	// Step-VI: Computation of Zernike moments starts 
	for(p=0; p<=pmax; p++)
	{
		pp1byarea=(p+1)/area;
		for(q=0; q<=p; q++)
		{
			if(!((p-q)%2==0)) continue;
			sumR=0.0;  sumI=0.0;  
            for(k=q; k<=p; k++)
			{
				if(!((p-k)%2==0)) continue;
				bpqk=B[p][q][k];
				sumR=sumR+bpqk*sR[q][k];
				sumI=sumI+bpqk*sI[q][k];
			} // k-loop
			sumR=pp1byarea*sumR;  sumI=pp1byarea*sumI;
            AnmReal[p][q]=sumR;
            AnmImg[p][q]=sumI;
			Anm[p][q]=sqrt(sumR*sumR+sumI*sumI);
		} // q-loop
	} // p-loop
	stop = clock(); // stop clock time
	fprintf(momptr, "\n    n   m    AnmReal    AnmImg     Anm        \n");
	for(p=0; p<=pmax; p++)
	{
		icount=0;
		for(q=0; q<=p; q++)
		{
			if(option==1)
			    if((p-abs(q))%2 !=0) continue;  // p-|q| must be even for Zernike
			if(q>=0)
			{
				if(icount==0)
				{
					fprintf(momptr, "%5d%5d%12.6f%12.6f%12.6f\n", 
						p, q, AnmReal[p][q], AnmImg[p][q], Anm[p][q]); 
					icount++;
				}
				else
					fprintf(momptr, "%10d%12.6f%12.6f%12.6f\n", 
					q, AnmReal[p][q], AnmImg[p][q], Anm[p][q]); 
			}  
		} // for q-loop
	}  // for p-loop
	cpu_time_used = ((double) (stop - start)) / CLOCKS_PER_SEC;
	fprintf(momptr, "Time taken by Zernike_moment_accurate() for the computation"
		" of moments is : %f sec\n", cpu_time_used); 
	return;
}

/*
/*************************************************************************
*                                                                        *
*   Calculation of Zernike moments through numerical integration. Zernike*
*   polynomials are calculated using q-recursive algorithm.	             *
*                                                                        *
**************************************************************************
*/
void Zernike_moment_numerical_integration(int N, int nmax, unsigned char fxy[MAXSIZE][MAXSIZE],
					double Anm[NMAX+1][2*NMAX+1], double AnmReal[NMAX+1][2*NMAX+1], 
					double AnmImg[NMAX+1][2*NMAX+1], FILE *momptr, unsigned char in_filename[14] )
{
	int n, m, p, q, x, y, icount, npoint, npixels=0, recursion=1;

	double mult, xf, yf, xfsqr, yfsqr, rsqr, rsqrt, 
		   pi=4.0*atan(1.0), Rn, Rnm, Rnm2, Rnmp2, Rnmp4, Rnnm2, H1, H2, H3,
		   np1byarea, a, b, C, S, temp1, temp2, fvalue, N2, wp, wq, tp, tq,  
		   rpower[NMAX+1], H1mn[NMAX+1][NMAX+1],  H2mn[NMAX+1][NMAX+1], 
		   H3mn[NMAX+1][NMAX+1], sinmt[NMAX+1], cosmt[NMAX+1],
		   w[MAXINT], t[MAXINT];
	N2=N*sqrt(2.0);  // outer circle
	printf("Computation of ZMs using numerical integration\n");
	fprintf(momptr, "\nOuter circle is used\n");
	Dfsqr=(N2/2.0)*(N2/2.0);
	start=clock(); // start clock time
	// printf("Enter no. of points for Gaussian integration : ");
	// scanf("%d", &npoint);
	// No. of point for Gaussian integration.
	// Maximum allowed value is 7.
	npoint = 5; 
	if(npoint>MAXINT) 
		display_error("No. of points for Gaussian integration more than allowed");
	if(npoint<1 || npoint>7) display_error("Wrong value for npoint for Gaussian integration");
	fprintf(momptr, "\n%dx%d-point Gaussian numerical integration is used\n", npoint, npoint);
	switch(npoint)
	{
		case 1 : w[0]=2.0;  t[0]=0.0;  break;
		case 2 : w[0]=1.0;  w[1]=w[0];  t[0]=-0.5773502692;  t[1]=-t[0];  break;
		case 3 : w[0]=0.5555555556;  w[1]=0.8888888889;   w[2]=w[0];
					t[0]=-0.7745966692;   t[1]=0.0;   t[2]=-t[0];   break;
		case 4 : w[0]=0.3478548451;   w[1]=0.6521451549;  w[2]=w[1];  w[3]=w[0];
					t[0]=-0.8611363116;  t[1]=-0.3399810436;   t[2]=-t[1];   t[3]=-t[0];  break;
		case 5 : w[0]=0.2369268851;   w[1]=0.4786286705;   w[2]=0.5688888889;   w[3]=w[1];   
					w[4]=w[0];
					t[0]=-0.9061798459;   t[1]=-0.5384693101;  t[2]=0.0;   t[3]=-t[1];   
					t[4]=-t[0];   break;
		case 6 : w[0]=0.1713244924;  w[1]=0.3607615730;  w[2]=0.4679139346;  w[3]=w[2];   
					w[4]=w[1];  w[5]=w[0];
					t[0]=-0.9324695142;  t[1]=-0.6612093865;  t[2]=-0.2386191861;  t[3]=-t[2];  
					t[4]=-t[1];  t[5]=-t[0];   break;
		case 7 : w[0]=0.1294849662;  w[1]=0.2797053915;   w[2]=0.3818300505;   w[3]=0.4179591837;
					w[4]=w[2];   w[5]=w[1];   w[6]=w[0];
					t[0]=-0.9491079123;   t[1]=-0.7415311856;   t[2]=-0.4058451514;   t[3]=0.0;   
					t[4]=-t[2];  t[5]=-t[1];   t[6]=-t[0];   break;
		default : display_error("Wrong value for the no. of points in Gaussian integration");
	}
	// Initialize the Zernike moments.
	for(n=0; n<=nmax; n++)
		for(m=0; m<=nmax; m++)
		{
			AnmReal[n][m]=0.0; AnmImg[n][m]=0.0;
			Anm[n][m]=0.0;
		}
	// compute H1, H2 and H3 and save. Values for m=n and m=n-2 are skipped
	for(n=0; n<=nmax; n++)
	{
		for(m=n-4; m>=0; m=m-2)
		{
			H3mn[n][m]=-(double)(4*(m+2)*(m+1))/((n+m+2)*(n-m));
			H2mn[n][m]=H3mn[n][m]*(n+m+4)*(n-m-2)/(4*(m+3))+(m+2);
			H1mn[n][m]=(m+4)*(m+3)/(double)2-H2mn[n][m]*(m+4)+
				H3mn[n][m]*(n+m+6)*(n-m-4)/8.0;
		}  // for m-loop
	}  // for n-loop
	npixels=0;
	for(p=0; p<npoint; p++)
	{
		wp=w[p];  tp=t[p];
		for(q=0; q<npoint; q++)
		{
			wq=w[q];  tq=t[q];
			for(x=0; x<N; x++)
			{
				xf=(tp+2.0*x+1-N)/N2;  xfsqr=xf*xf;
				for(y=0; y<N; y++)
				{
					yf=(tq+2.0*y+1-N)/N2;  yfsqr=yf*yf;
					// skip those integration sampling points which fall outside the 
					// circle to make sure that |r|<=1, otherwise value of radial 
					// polynomials will be very high. This condition is usefule when 
					// we remove numerical integration error and geometric error without
					// using arc-grids as the approximation of the circular boundary will
					// be better.
					rsqrt=sqrt(xfsqr+yfsqr);
					if(rsqrt>1.0) continue;
					npixels++;
					fvalue=fxy[x][y];
					rsqr=rsqrt*rsqrt;
					// Computue r**n
					mult=rsqrt;
					rpower[0]=1.0;
					for(n=1; n<=nmax; n++)
					{
						rpower[n]=mult;
						mult=mult*rsqrt;
					}  // for n-loop
					cosmt[0]=1.0;  	sinmt[0]=0.0;
					// for m=1
					if(rsqrt==0.0) { a=1.0;  b=0.0; }  // theta is 0.0 for xf=0.0 & yf=0.0
					else { a=xf/rsqrt;  b=yf/rsqrt; }
					cosmt[1]=a;  
					sinmt[1]=b;
					C=a;  S=b;
					// for m>=2
					for(m=2; m<=nmax; m++)
					{
						temp1=a*C-b*S;
						temp2=a*S+b*C;
						cosmt[m]=temp1;  
						sinmt[m]=temp2;
						C=temp1;  S=temp2;
					}  // m-loop
					for(n=0; n<=nmax; n++)
					{
						Rn=rpower[n];
						if(n>1) Rnm2=rpower[n-2];
						for(m=n; m>=0; m=m-2)
						{
							if(m==n)
							{
								Rnm=Rn;
								Rnmp4=Rn;  // save Rn(m+4)
							}
							else if(m==n-2)
							{
								Rnnm2=n*Rn-(n-1)*Rnm2;
								Rnm=Rnnm2;
								Rnmp2=Rnnm2;  // save Rn(m+2)
							}
							else
							{
								H3=H3mn[n][m];
								H2=H2mn[n][m];
								H1=H1mn[n][m];
								if(rsqr==0.0)
								{
									if(m==0)
									{
										if(n%4==0) Rnm=1.0;
										else Rnm=-1.0;
									}
									else Rnm=0.0;
								}  // if-rsqr
								else Rnm=H1*Rnmp4+(H2+H3/rsqr)*Rnmp2;
								Rnmp4=Rnmp2;  // Rnmp2 now becomes Rnmp4 for next m
								Rnmp2=Rnm;  // Rnm becomes Rnmp2 for next m
							}  // if block
							// printf("n=%d, m=%d, Rnm=%f\n", n, m, Rnm);  getch();
							AnmReal[n][m]=AnmReal[n][m]+0.25*wp*wq*fvalue*Rnm*cosmt[m];
							AnmImg[n][m]=AnmImg[n][m]+0.25*wp*wq*fvalue*Rnm*sinmt[m];
						}  // for m-loop
					}  // for n-loop
				}  // for y-loop
			}  // for x-loop
		}  // for q-loop
	}  // for p-loop
	stop = clock(); // stop clock time
	fprintf(momptr, "\n    n   m    AnmReal    AnmImg     Anm        \n");
	double area=pi*Dfsqr;
	for(n=0; n<=nmax; n++)
	{
		np1byarea=(n+1)/area;
		icount=0;
		for(m=0; m<=n; m++)
		{
			if((n-abs(m))%2 !=0) continue;  // n-|m| must be even for Zernike
			AnmReal[n][m]=np1byarea*AnmReal[n][m];
			AnmImg[n][m]=-np1byarea*AnmImg[n][m];
			Anm[n][m] = sqrt(AnmReal[n][m]*AnmReal[n][m]+AnmImg[n][m]*AnmImg[n][m]);
			if(m>=0)
			{
				if(icount==0)
				{
					fprintf(momptr, "%5d%5d%12.6f%12.6f%12.6f\n", 
						n, m, AnmReal[n][m], AnmImg[n][m], Anm[n][m]); 
					icount++;
				}
				else
					fprintf(momptr, "%10d%12.6f%12.6f%12.6f\n", 
					m, AnmReal[n][m], AnmImg[n][m], Anm[n][m]); 
			}  
		} // for m-loop
	}  // for n-loop
	cpu_time_used = ((double) (stop - start)) / CLOCKS_PER_SEC;
	fprintf(momptr, "Time taken by Zernike_moment_numerical_integration() for"
		" the computation of moments is : %f sec\n", cpu_time_used);
	return;
}

/*
/*************************************************************************
*                                                                        *
*   Reconstruction of the image using computed ZMs.                      * 
*                                                                        *
**************************************************************************
*/
void Reconstruction_ZM(int N, int nmax, double Anm[NMAX+1][2*NMAX+1], 
	double AnmReal[NMAX+1][2*NMAX+1], double AnmImg[NMAX+1][2*NMAX+1], 
	unsigned char fxyoriginal[MAXSIZE][MAXSIZE],
	unsigned char fxyout[MAXSIZE][MAXSIZE], FILE *momptr, unsigned char in_filename[14])
{
	int n, m, x, y, ndiff=0, nblack, npixels=0;
	double mult, xf, yf, xfsqr, yfsqr, rsqr, rsqrt, pi=4.0*atan(1.0), sum,   
		   Rn, Rnm, Rnm2, Rnmp2, Rnmp4, Rnnm2, H1, H2, H3, funvalue, cost, sint, 
		   a, b, C, S, temp, N2, rpower[NMAX+1], H1mn[NMAX+1][NMAX+1],  H2mn[NMAX+1][NMAX+1], 
		   H3mn[NMAX+1][NMAX+1], fxy_double[MAXSIZE][MAXSIZE], fxyR[MAXSIZE][MAXSIZE], 
		   cosmt[NMAX+1], sinmt[NMAX+1]; 
	N2=N*sqrt(2.0);  // outer circle
	Dfsqr=(N2/2.0)*(N2/2.0);
	// Initialize the functions.
	for(x=0; x<MAXSIZE; x++)
	{
		for(y=0; y<MAXSIZE; y++)
		{
			fxyR[x][y]=255; fxyout[x][y]=255;
		}
	}
	// Copy image data to real values
	for(x=0; x<N; x++)
	{
		for(y=0; y<N; y++)
		{
			fxy_double[x][y]=fxyoriginal[x][y];
		}
	}
	// printf("Calling the function-Reconstrion");  getch();
	// Compute H1, H2 and H3 for q-recursive and save. Values for m=n and m=n-2 are skipped
	for(n=0; n<=nmax; n++)
	{
		for(m=n-4; m>=0; m=m-2)
		{
			H3mn[n][m]=-(double)(4*(m+2)*(m+1))/((n+m+2)*(n-m));
			H2mn[n][m]=H3mn[n][m]*(n+m+4)*(n-m-2)/(4*(m+3))+(m+2);
			H1mn[n][m]=(m+4)*(m+3)/(double)2-H2mn[n][m]*(m+4)+
				H3mn[n][m]*(n+m+6)*(n-m-4)/8.0;
		}  // for m-loop
	}  // for n-loop
	for(x=0; x<N; x++)
	{
		xf=(2.0*x+1-N)/N2;  xfsqr=xf*xf;
		for(y=0; y<N; y++)
		{
			yf=(2.0*y+1-N)/N2;  yfsqr=yf*yf;
			rsqrt=sqrt(xfsqr+yfsqr);
			if(rsqrt>1.0) continue;
			npixels++;
			rsqr=rsqrt*rsqrt;
			// Computue r**n
			mult=rsqrt;
			rpower[0]=1.0;  
			for(n=1; n<=nmax; n++)
			{
				rpower[n]=mult;
				mult=mult*rsqrt;
			}  // for n-loop
			// Comput cos(m*theta) and sin(m*theta)
			// For m=0 and m=1
			cosmt[0]=1.0;  	sinmt[0]=0.0;
			if(rsqrt==0.0) { a=1.0;  b=0.0; }  // theta is 0.0 for xf=0.0 & yf=0.0
			else { a=xf/rsqrt;  b=yf/rsqrt; }
			cosmt[1]=a;  
			sinmt[1]=b;
			C=a;  S=b;
			// For m>=2
			for(m=2; m<=nmax; m++)
			{
				temp=a*C-b*S;
				S=a*S+b*C;
				cosmt[m]=temp;  
				sinmt[m]=S;
				C=temp;  
			}  // for m-loop
			funvalue=0.0;
			for(n=0; n<=nmax; n++)
			{
				Rn=rpower[n];
				if(n>1) Rnm2=rpower[n-2];
				for(m=n; m>=0; m=m-2)
				{
					if(m==n)
					{
						Rnm=Rn;
						Rnmp4=Rn;  // save Rn(m+4)
					}
					else if(m==n-2)
					{
						Rnnm2=n*Rn-(n-1)*Rnm2;
						Rnm=Rnnm2;
						Rnmp2=Rnnm2;  // save Rn(m+2)
					}
					else
					{
						H3=H3mn[n][m];
						H2=H2mn[n][m];
						H1=H1mn[n][m];
						if(rsqr==0.0)
						{
							if(m==0)
							{
								if(n%4==0) Rnm=1.0;
								else Rnm=-1.0;
							}
							else Rnm=0.0;
						}  // if-rsqr
						else Rnm=H1*Rnmp4+(H2+H3/rsqr)*Rnmp2;
						Rnmp4=Rnmp2;  // Rnmp2 now becomes Rnmp4 for next m
						Rnmp2=Rnm;  // Rnm becomes Rnmp2 for next m
					}  // if block
					cost=cosmt[m];
					sint=sinmt[m];
					if(m==0) sum = Rnm*(AnmReal[n][m]*cost - AnmImg[n][m]*sint);
					else sum=2.0*Rnm*(AnmReal[n][m]*cost	- AnmImg[n][m]*sint);
					// printf("n=%d, m=%d, Anm=%f, sum=%f\n", n, m, Anm[n][m], sum);  
					// getch();
					funvalue = funvalue+sum;  
				}  // for m-loop
			}  // for n-loop
			// printf("x=%d, y=%d, N2=%d, rsqr=%f, funvalue=%f, fxyoriginal[x][y]=%d\n", 
			//	x, y, N2, rsqr, funvalue, fxyoriginal[x][y]);  getch(); 
			fxyR[x][y]=funvalue;
		}  // for y-loop
	}  // for x-loop
	npixels=0;  nblack=0;
	for(x=0; x<N; x++)
	{
		xf=(2.0*x+1-N)/N2;
		for(y=0; y<N; y++)
		{
		    yf=(2.0*y+1-N)/N2;
			if((xf*xf+yf*yf)>1.0) continue;
			npixels++;
			if(fxyR[x][y]<0.0) fxyout[x][y]=0;
			else if(fxyR[x][y]>255.0) fxyout[x][y]=255;
			else fxyout[x][y]=(unsigned char)fxyR[x][y];
		}  // for y-loop
	}  // for x-loop
	MSRE(N, fxyoriginal, fxyR, momptr);
}

/* Calculates Mean square reconstruction error(MSRE)  */
/*************************************************************************
*                                                                        *
*   Calculation of Mean square reconsturction error.                     *
*                                                                        *
**************************************************************************
*/
void MSRE(int N, unsigned char fxyoriginal[MAXSIZE][MAXSIZE], 
		  double fxyR[MAXSIZE][MAXSIZE], FILE *momptr)
{
	int x, y, npixels;
	unsigned char fxyij;
	double N2, msre=0.0, rho=0.0, xf, yf, xfsqr, yfsqr, rsqr, sum1, sum2, diff, fxycij;
	N2=N*sqrt(2.0);  // outer circle, for inner circle use N2 = N
	npixels=0;  sum1=0; sum2=0;
	for(x=0; x<N; x++)
	{
		xf=(2*x+1-N)/N2;  xfsqr=xf*xf;
		for(y=0; y<N; y++)
		{
            yf=(2*y+1-N)/N2;   yfsqr=yf*yf;
			rsqr=xfsqr+yfsqr;
			// Use outer circle for MSRE calculation
			if(rsqr>1.0) continue;  
			fxyij=fxyoriginal[x][y];
			fxycij=fxyR[x][y];
			diff=(double)fxyij-fxycij;
			sum1=sum1+(double)fxyij*(double)fxyij;
			sum2=sum2+diff*diff;
		}  // y-loop
	}  // x-loop
	if(sum2==0.0) fprintf(momptr, "Output image is the same as the input image\n");
	if(sum1==0.0) display_error("Original image is blank");
	else msre=sum2/sum1;
	printf("Mean Square Reconstruction ZM Error(eps) : %f\n", msre);  
	fprintf(momptr, "Mean Square Reconstruction ZM Error(eps) : %10.5f", msre);
}

// Copies the gray scale image to a 2D array
/*************************************************************************
*                                                                        *
*  Copies the grayscale image to a 2D array.                             *
*                                                                        *
**************************************************************************
*/
void copy_to_array(unsigned long ImageLength, unsigned long ImageWidth, 
				   unsigned char *IB, unsigned char fxy[MAXSIZE][MAXSIZE])
{
	   int npixels=0; 
	   unsigned int i, j;
	   for(i=0; i<ImageLength; i++)
	   {
		   for(j=0; j<ImageWidth; j++)
		   {
			   npixels++;
			   fxy[i][j]=*IB;
			   IB++;
		   }  //  j-loop
	   }  // i-loop
}

// Copies the array to image output file
/*************************************************************************
*                                                                        *
*  Copies the array to image output file                                 *
*                                                                        *
**************************************************************************
*/
void copy_array_to_file(unsigned long ImageLength, unsigned long ImageWidth, 
				   unsigned char fxyout[MAXSIZE][MAXSIZE], FILE *outfile_ptr)
{
	   int npixels=0; 
	   unsigned int i, j;
	   unsigned char pic_byte;
	   for(i=0; i<ImageLength; i++)
	   {
		   for(j=0; j<ImageWidth; j++)
		   {
			   npixels++; pic_byte=fxyout[i][j];  // R-component
			   fwrite((const void *) &pic_byte, (size_t)sizeof(pic_byte), (size_t)1, 
				   outfile_ptr);
		   }  //  j-loop
	   }  // i-loop
	   // printf("For Output Image : Total pixels=%d\n", npixels);  getch();   
}

// Displays the error messages.
/*************************************************************************
*                                                                        *
*   Displays error.                                                      *
*                                                                        *
**************************************************************************
*/
void display_error( char *message)
{
	printf("%s\n", message);
	printf("Program stopped !!\n");  getch();
	exit(0);
}