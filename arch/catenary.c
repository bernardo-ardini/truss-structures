double lcat;

void Tcat(struct matrix *a,struct matrix *t)
{
	loadm(t,1,1);
	double A=getm(a,0,0);
	setm(t,0,0,A*(cosh(lcat/A)-1));
}

void Kcat(struct matrix *a,struct matrix *k)
{
	loadm(k,1,1);
	double A=getm(a,0,0);
	setm(k,0,0,cosh(lcat/A)-1-lcat/A*sinh(lcat/A));
}

double find_a(double h,double l)
{
	lcat=l;
	
	struct matrix A;
	zeros(&A,1,1);
	setm(&A,0,0,lcat/asinh(2*h/l));
	
	struct matrix H;
	loadm(&H,1,1);
	
	setm(&H,0,0,h);
	
	struct matrix a;
	
	nr(&A,&a,Kcat,Tcat,&H,1e-6);
	
	double res=getm(&a,0,0);
	
	clearm(&A);
	clearm(&H);
	clearm(&a);
	
	return res;
}