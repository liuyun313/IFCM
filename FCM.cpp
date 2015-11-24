#include "iostream"
#include "FCM.h"
#include "memory.h"
#include "stdlib.h"
#include "stdio.h"
#include "cstring"
#include "math.h"
#include "time.h"

namespace ker
{
	#define CONST_E 1e-5

	FCM::FCM(int a,int b,int c)
	{
		nVec=a;
		nDim=b;
		nClu=c;
		J=0;
		q=2;

		//membership matrix initialization using random number 
		u=new double*[nVec];
		for(int i=0;i<nVec;i++)
		{
			u[i]=new double[nClu];
			int n[nClu],sum=0;
			//srand((unsigned)time(NULL)); 
			for(int j=0;j<nClu;j++)
			{
				n[j]=rand()%100;
				sum+=n[j];
			}
			for(int j=0;j<nClu;j++)u[i][j]=(double)n[j]/sum;			
		}

		//vLabel initialization
		vLabel=new int[nVec];
		memset(vLabel,0,sizeof(int)*nVec);
		int lab;
		for(int i=0;i<nVec;i++)
		{
			double m=0;
			for(int j=0;j<nClu;j++)
			{
				if(u[i][j]>m)
				{
					vLabel[i]=j;
					m=u[i][j];
				}
			}
		}
//for(int i=0;i<nClu;i++)for(int j=0;j<nVec;j++)printf("%.5f ",pow(u[j][i],q));

		//cluster-size initialization
		cSize=new double[nVec];
		memset(cSize,0,sizeof(double)*nVec);		
		for(int i=0;i<nClu;i++)
		{
			for(int j=0;j<nVec;j++)cSize[i]+=u[j][i];
			cSize[i]/=nVec;
		}

		//clustering center initialization
		center=new double*[nClu];
		for(int i=0;i<nClu;i++)
		{
			center[i]=new double[nDim];
			memset(center[i],0,sizeof(double)*nDim);
		}
	}

	FCM::~FCM()
	{
		for(int i=0;i<nVec;i++)delete u[i];
		delete []u;
		for(int i=0;i<nClu;i++)delete center[i];
		delete []center;
	}
	
	double FCM::dist(double *a,double *b)
	{
		double d=0;
		for(int i=0;i<nDim;i++)
		{
			d+=pow((*(a+i)-*(b+i)),2);
		}
		d=sqrt(d);
		return d;
	}

	void FCM::updateCenter()//update clustering center using membership matrix
	{

		double a=0,b=0;
//for(int k=0;k<nVec;k++)for(int i=0;i<nClu;i++){b+=pow(u[k][i],q);printf("%.5f ",b);}
		for(int j=0;j<nDim;j++)
		{
			for(int i=0;i<nClu;i++)
			{
				a=0;
				b=0;
				for(int k=0;k<nVec;k++)
				{
					a+=pow(u[k][i],q)*fVector[k][j];
					b+=pow(u[k][i],q);
				}
				center[i][j]=a/b;//
			}
			
		}
	}

	
	void FCM::updateU()//update membership matrix 
	{
		double d1[nClu],d2=0,u1[nVec][nClu];
		for(int i=0;i<nVec;i++)
			for(int j=0;j<nClu;j++)u1[i][j]=u[i][j];


		for(int i=0;i<nVec;i++)
		{
			d2=0;
			for(int j=0;j<nClu;j++)
			{
				d1[j]=dist(fVector[i],center[j]);
				d1[j]=pow(d1[j],2/(q-1));
				if(d1[j]<1e-20)d1[j]=1e-20;
				d2+=1/d1[j];
			}
			for(int j=0;j<nClu;j++)
			{
				u[i][j]=cSize[i]/d1[j]/d2;
			}
		}


		//cSize calculation 
		//cSize=new double[nVec];
		memset(cSize,0,sizeof(double)*nVec);
		for(int i=0;i<nClu;i++)
		{
			for(int j=0;j<nVec;j++)cSize[i]+=u[j][i];
			cSize[i]/=nVec;
		}


		//calculate price value J
		J=0;
		for(int i=0;i<nVec;i++)
		{
			for(int j=0;j<nClu;j++)
			{
				if(fabs(u[i][j]-u1[i][j])>J)J=fabs(u[i][j]-u1[i][j]);
			}
		}
	}


	void FCM::Cluster(double **X,int *label)
	{


		fVector=X;

		double J1;
		do
		{
			J1=J;
			updateCenter();
			updateU();
		}while(fabs(J-J1)>CONST_E);
//for(int i=0;i<nDim;i++)printf("%.5f ",center[0][i]);
//for(int i=0;i<nClu;i++)printf("%.5f ",u[0][i]);
		FS=0;
		double *center_;
		center_=new double [nDim];
		memset(center_,0,sizeof(double)*nDim);
		for(int i=0;i<nDim;i++)
		{
			for(int j=0;j<nClu;j++)center_[i]+=center[j][i];
			center_[i]/=nClu;
		}

		for(int j=0;j<nClu;j++)
		{
			FS-=dist(center[j],center_)/nClu;
			for(int i=0;i<nVec;i++)
			{
				FS+=pow(u[i][j],q)*dist(fVector[i],center[j])/nClu/nVec;
			}
		}
		//FS/=nClu;
//printf("%.5f\n",dist(fVector[0],fVector[1]));
		for(int i=0;i<nVec;i++)*(label+i)=vLabel[i];






























		

		
	}









}

