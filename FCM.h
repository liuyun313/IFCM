#ifndef FCM_H
#define FCM_H

namespace ker
{
	class FCM
	{
	private:
		int nVec;
		int nDim;
		int nClu;
		double q;//fuzzy degree
		double FS;//modified PC

		double **fVector;//feature vector
		double **u;//membership matrix
		double **center;//clustering centers
		double *cSize;
		int *vLabel;//cluster label of each point
		double J;//price value	
	public:
		FCM(int a,int b,int c);
		~FCM();
		void Cluster(double **X,int *label);
		double getJ(){return J;}
		double getF(){return FS;}
	private:
		void updateU();
		void updateCenter();
		double dist(double *a,double *b);
	};
}
#endif
