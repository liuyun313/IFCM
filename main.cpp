#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<cstring>
#include<fstream>
#include<memory.h>
#include<math.h>
#include "FCM.h"
//#include"class.h"
#include<time.h>
//#include"FCM.cpp"



using namespace ker;
using namespace std;

//******************************end****************************//

void kmer(char *s,int *p,int bp)//kmer calculation****n reads with l bp
{
    int k=4,l=bp;
    int *g = new int[l];
    memset(g,0,sizeof(int)*l);
    int fVector[256]={0};
	//********convert the "ATGC" to "0123"
	for(int i=0;i<l;i++)
	{
		if(*s=='A')g[i]=0;
		if(*s=='T')g[i]=1;
		if(*s=='G')g[i]=2;
		if(*s=='C')g[i]=3;
		s++;
	}
	//*****kmer calculate*****
	for(int i=0;i<l-3;i++)
	{
		fVector[64*g[i]+16*g[i+1]+4*g[i+2]+g[i+3]]++;		
	}
	delete g;
	//***********k-mer reduction
	int weight[4]={1,4,16,64};
	for(int i=0;i<256;i++)
	{
		if(fVector[i]==-1)
			continue;
		else		
		{
			int a,b,c,d;
			a=(int)(i/64);
			b=(int)((i-a*64)/16);
			c=(int)((i-a*64-b*16)/4);
			d=i-a*64-b*16-c*4;
			int rank[4]={a,b,c,d};
			for(int j=0;j<4;j++)
			{
				if(rank[i]==0){rank[i]=1;continue;}
				if(rank[i]==1){rank[i]=0;continue;}
				if(rank[i]==2){rank[i]=3;continue;}
				if(rank[i]==3){rank[i]=2;continue;}
			}
			int i1=0;
			for(int j=0;j<4;j++)
				i1+=rank[j]*weight[j];
			if(i1!=i)
			{
				fVector[i]+=fVector[i1];
				fVector[i1]=-1;
			}
		}
	}

	int vec[136]={0},num=0;
	
	for(int i=0;i<256;i++)
	{
		if(fVector[i]!=-1)
		{
			vec[num]=fVector[i];
			num++;
		}
	}

	//return fVector;
	for(int i=0;i<136;i++)
	{	
		*(p+i)=vec[i];
	}


}





void normalization(int **X1,double **X,int num)
{
	double max,min;
	for(int i=0;i<136;i++)
	{
		max=X1[0][i];
		for(int j=1;j<num;j++)
		{
			if(X1[j][i]>max)max=X1[j][i];
		}
		for(int j=0;j<num;j++)
		{
			if(max!=0)X[j][i]=X1[j][i]/max;
		}
	}
}



int main(int argc, char *argv[])
{
	time_t t1,t2;
	t1=time(NULL);
	char f[100];
	strcpy(f, argv[1]);
	
	
//***********************file read*******************//
	//calculate the number of symbols in the input file
	ifstream in(f);//("1a.fna");
	in.seekg(0, ios::end);      //set the file point to the end of the file stream
	streampos ps = in.tellg();  //read the position the poing,while ps will be the symbol number
	in.close();                 //鍏抽棴鏂囦欢娴?
	FILE *fd;
	int len=int(ps);
	char *str=new char[len];
	fd=fopen(f,"r");//("1a.fna","r");
	fread(str,1,len,fd);//read the file into str
	fclose(fd);
//******************************end*******************//

//************separate the title and the DNA read according to '\n'********//
	int num=0;
	for(int i=0;i<len;i++)
	{
		if(str[i]=='>')num++;
	}

	int bp[num];//length of each read
	num=0;
	int i=0;
	while(i<len)
	{
		if(str[i]=='>')
		{
			int j=0;
			while(!(str[i]=='\n'))
			{
				i++;
				j++;
			}
			i++;
			j=0;
			while(!(str[i]=='>')&&i<len)
			{
			   if(!(str[i]=='\n'))
			   {
				j++;
			   }	
				i++;
			}
			bp[num]=j-1;
			i--;
			num++;
		}
		i++;
	} 

	int length=0;
	for(i=0;i<num;i++)
		if(bp[i]>length)length=bp[i];

       char **title=new char*[num];
	for(int i=0;i<num;i++)
	{
		title[i]=new char[100]; 
		memset(title[i],0,sizeof(char)*100);
	}
	
       char **read= new char*[num];
	for(int i=0;i<num;i++)
	{
		read[i]=new char[length];
		memset(read[i],0,sizeof(char)*length);
	}

///char title[num][length]=new char*[num*length] ,read[num][length];

	num=0;	
	while(i<len)
	{
		if(str[i]=='>')
		{
			int j=0;
			while(!(str[i]=='\n'))
			{
				title[num][j]=str[i];
				i++;
				j++;
			}
			i++;
			j=0;
			while(!(str[i]=='>')&&i<len)
			{
			   if(!(str[i]=='\n'))
			   {
			       read[num][j]=str[i];
				j++;
			   }	
				i++;
			}
			bp[num]=j-1;
			i--;
			num++;
		}
		i++;
	} 

	delete str;
//*********************end**********************//
        
	cout<<'a'<<endl;
//***************calculate the k-mer frequency*************//
	//int x[256]={0},*p;
	int **X1=new int*[num];//feature vector with dimension num*256
	for(i=0;i<num;i++)
	{
		X1[i]=new int[136];
		memset(X1[i],0,sizeof(int)*136);
	}
	//p=x;
	for(i=0;i<num;i++)
	{
		kmer(read[i],X1[i],bp[i]);		
	}
//********************end******************//


//**************feature vector normalization******//
	

	double **X = new double*[num];
	for (int i = 0; i < num; i++)
	{
		X[i] = new double[136];
	}

	normalization(X1,X,num);
	
//////////////////////////////////////////////////////
	for(i=0;i<num;i++)delete X1[i];
	delete[] X1;
//for(i=0;i<256;i++)cout<<X[0][i]<<' ';
//******************end*********************//

//*****************species number estimation****************//
	int l=0,kmin=0,kmax=0;//sequencing depth, average length of DNA reads and estimated species number
	int Gmin=1.18e6,Gmax=5.24e6;
	for(i=0;i<num;i++)l+=bp[i];
	kmin=(int)l/Gmax;
	kmax=(int)l/Gmin;
	if(kmin<3){kmin=3;kmax=5*kmin;}
	cout<<"The number of contigs: "<<num<<endl;
	cout<<"The total lenth of contigs :"<<l<<endl;
	cout<<"kmin="<<kmin<<endl;
	cout<<"kmax="<<kmax<<endl;

//********************end********************//


//******c-means clustering for the initials of model-based method****//
	int n=kmax-kmin+1;
	int **label;//[n][num];
	label=new int*[n];
	for(i=0;i<n;i++)
	{
		label[i]=new int[num];
		memset(label[i],0,sizeof(int)*num);
	}

	//cmeans(p_X,p_label,num,kmax,kmin);
//************************end**************//

//*****************FCM clustering method********************//
	
	double J[n],Jmin;
	for(i=0;i<n;i++)
	{
		FCM *fcm=new FCM(num,136,kmin+i);
		fcm->Cluster(X,label[i]);//cout<<X[0][0]<<endl;
		J[i]=fcm->getF();
		cout<<kmin+i<<' '<<J[i]<<endl;
		delete fcm;
	}

	Jmin=J[0];
	int k_optimal=kmin;
	for(i=1;i<n;i++)
	{
		if(J[i]<Jmin)
		{
			Jmin=J[i];
			k_optimal=i+kmin;
		}
	}
	//cout<<"The number of optimal number of clusters is "<<k_optimal<<endl;

//*******************end********************//	
/*FCM fcm(num,256,3);
fcm.Cluster(X,label[0]);*/

	//int np[k_optimal];
	//for(i=0;i<k_optimal;i++)np[i]=0;
	//for(i=0;i<num;i++)np[label[k_optimal-kmin][i]]++;
	//for(i=0;i<k_optimal;i++)cout<<np[i]<<endl;



//output the binning result
	cout<<"Writing binning result into file......"<<endl;
	ofstream out("result.txt");
	if(out.is_open())
	{
		//out<<"the cluster label of each read\n";
		for(int j=0;j<k_optimal;j++)
		{
			out<<"Clusters "<<j<<endl;
			for(i=0;i<num;i++)
			{
				if(label[k_optimal-kmin][i]==j)out<<title[i]<<endl;								
			}
		}
		out.close();
	}

	for(int i=0;i<num;i++)
	{
		delete[] title[i];
		delete[] read[i];
	}
	delete[] title;
	delete[] read;
 
	t2=time(NULL);
	cout<<"The running time is "<<difftime(t2,t1)<<'s'<<endl;



//****************end*****************//
}


