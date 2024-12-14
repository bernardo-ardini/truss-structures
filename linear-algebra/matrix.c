struct matrix
{
	double *pointer;
	int n;
	int m;
};

void loadm(struct matrix *A,int n,int m);
void clearm(struct matrix *A);
void row(struct matrix *A,int i,struct matrix *r);
void column(struct matrix *A,int j,struct matrix *c);
double norm(struct matrix *v);
void difference(struct matrix *A,struct matrix *B,struct matrix *C);
void _difference(struct matrix *A,struct matrix *B);
void sum(struct matrix *A,struct matrix *B,struct matrix *C);
void _sum(struct matrix *A,struct matrix *B);
void product(double a,struct matrix *A,struct matrix *B);
void _product(double a,struct matrix *A);
void setm(struct matrix *A,int i,int j,double a);
double getm(struct matrix *A,int i,int j);
void printm(struct matrix *A);
void matmul(struct matrix *A,struct matrix *B,struct matrix *C);
void dyadic(struct matrix *v,struct matrix *vv);
void eye(struct matrix *A,int n);
void zeros(struct matrix *A,int n,int m);
void copy(struct matrix *A,struct matrix *B);

void loadm(struct matrix *A,int n,int m)
{
	A->n=n;
	A->m=m;
	A->pointer=(double*)malloc(n*m*sizeof(double));
}

void clearm(struct matrix *A)
{
	A->n=0;
	A->m=0;
	free(A->pointer);
	A->pointer=NULL;
}

void row(struct matrix *A,int i,struct matrix *r)
{
	loadm(r,1,A->m);
	for(int j=0;j<(A->m);j++)
	{
		*((r->pointer)+j)=*((A->pointer)+i*A->m+j);
	}
}

void column(struct matrix *A,int j,struct matrix *c)
{
	loadm(c,A->n,1);
	for(int i=0;i<(A->n);i++)
	{
		*((c->pointer)+i)=*((A->pointer)+i*(A->m)+j);
	}
}

double norm(struct matrix *v)
{
	double sum=0;
	for(int i=0;i<(v->n);i++)
	{
		sum+=pow(*((v->pointer)+i),2);
	}
	return sqrt(sum);
}

void difference(struct matrix *A,struct matrix *B,struct matrix *C)
{
	loadm(C,A->n,A->m);
	for(int i=0;i<(A->n);i++)
	{
		for(int j=0;j<(A->m);j++)
		{
			*(C->pointer+i*A->m+j)=*(A->pointer+i*A->m+j)-*(B->pointer+i*A->m+j);
		}
	}
}

void _difference(struct matrix *A,struct matrix *B)
{
	for(int i=0;i<(A->n);i++)
	{
		for(int j=0;j<(A->m);j++)
		{
			*(A->pointer+i*A->m+j)-=*(B->pointer+i*A->m+j);
		}
	}
}

void sum(struct matrix *A,struct matrix *B,struct matrix *C)
{
	loadm(C,A->n,A->m);
	for(int i=0;i<(A->n);i++)
	{
		for(int j=0;j<(A->m);j++)
		{
			*(C->pointer+i*A->m+j)=*(A->pointer+i*A->m+j)+*(B->pointer+i*A->m+j);
		}
	}
}

void _sum(struct matrix *A,struct matrix *B)
{
	for(int i=0;i<(A->n);i++)
	{
		for(int j=0;j<(A->m);j++)
		{
			*(A->pointer+i*A->m+j)+=*(B->pointer+i*A->m+j);
		}
	}
}

void product(double a,struct matrix *A,struct matrix *B)
{
	loadm(B,A->n,A->m);
	for(int i=0;i<(A->n);i++)
	{
		for(int j=0;j<(A->m);j++)
		{
			*(B->pointer+i*A->m+j)=a**(A->pointer+i*A->m+j);
		}
	}
}

void _product(double a,struct matrix *A)
{
	for(int i=0;i<(A->n);i++)
	{
		for(int j=0;j<(A->m);j++)
		{
			*(A->pointer+i*A->m+j)*=a;
		}
	}
}

void setm(struct matrix *A,int i,int j,double a)
{
	*((A->pointer)+(A->m)*i+j)=a;
}

double getm(struct matrix *A,int i,int j)
{
	return *((A->pointer)+(A->m)*i+j);
}

void printm(struct matrix *A)
{
	for(int i=0;i<(A->n);i++)
	{
		for(int j=0;j<(A->m);j++)
		{
			//if(j>0 && j%(A->m)==0) printf("\n");
			printf("%.2f\t",getm(A,i,j));
		}
		printf("\n");
	}
}

void fprintm(char *filename,struct matrix *A,char *mode)
{
	FILE *file;
	file=fopen(filename,mode);
	
	for(int i=0;i<(A->n);i++)
	{
		for(int j=0;j<(A->m);j++)
		{
			if(j>0 && j%(A->m)==0) printf("\n");
			fprintf(file,"%.2f\t",getm(A,i,j));
		}
		fprintf(file,"\n");
	}
	
	fclose(file);
}

void matmul(struct matrix *A,struct matrix *B,struct matrix *C)
{
	int n=A->n;
	int m=A->m;
	int l=B->m;
	loadm(C,n,l);
	
	for(int i=0;i<n*l;i++)
	{
		*(C->pointer+i)=0;
	}

	#pragma omp parallel for
	for(int i=0;i<n;i++)
	{
		for(int k=0;k<m;k++)
		{
			for(int j=0;j<l;j++)
			{
				*(C->pointer+i*l+j)+=*(A->pointer+i*m+k)**(B->pointer+k*l+j);
			}
		}
	}
}

void dyadic(struct matrix *v,struct matrix *vv)
{
	int n=v->n;
	loadm(vv,n,n);
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			*(vv->pointer+i*n+j)=*(v->pointer+i)**(v->pointer+j);
		}
	}
}

void eye(struct matrix *A,int n)
{
	loadm(A,n,n);
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			if(i==j) *(A->pointer+i*n+j)=1;
			else *(A->pointer+i*n+j)=0;
		}
	}
}

void zeros(struct matrix *A,int n,int m)
{
	loadm(A,n,m);
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
			*(A->pointer+i*m+j)=0;
		}
	}
}

void ones(struct matrix *A,int n,int m)
{
	loadm(A,n,m);
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
			*(A->pointer+i*m+j)=1;
		}
	}
}

void copy(struct matrix *A,struct matrix *B)
{
	int n=A->n;
	int m=A->m;
	loadm(B,n,m);
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
				*(B->pointer+i*m+j)=*(A->pointer+i*m+j);
		}
	}
}

void _append_column(struct matrix *A,struct matrix *c)
{
	int n=A->n;
	int m=A->m;
	
	struct matrix B;
	loadm(&B,n,m+1);
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
				*(B.pointer+i*(m+1)+j)=*(A->pointer+i*m+j);
		}
	}
	
	clearm(A);
	
	for(int i=0;i<n;i++)
	{
		*(B.pointer+i*(m+1)+m)=*(c->pointer+i);
	}

	A->pointer=B.pointer;
	A->n=n;
	A->m=m+1;
}

void _append_row(struct matrix *A,struct matrix *r)
{
	int n=A->n;
	int m=A->m;
	
	struct matrix B;
	loadm(&B,n+1,m);
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
				*(B.pointer+i*m+j)=*(A->pointer+i*m+j);
		}
	}
	
	clearm(A);
	
	for(int j=0;j<m;j++)
	{
		*(B.pointer+n*m+j)=*(r->pointer+j);
	}

	A->pointer=B.pointer;
	A->n=n+1;
	A->m=m;
}

void readsubm(struct matrix *A,struct matrix *B,int i1,int i2,int j1,int j2)
{
	int n=A->n;
	int m=A->m;

	int N=i2-i1+1;
	int M=j2-j1+1;
	
	loadm(B,N,M);
	
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			*(B->pointer+M*i+j)=*(A->pointer+m*(i+i1)+j+j1);
		}
	}
}

void _writesubm(struct matrix *A,struct matrix *B,int i1,int i2,int j1,int j2)
{
	int n=A->n;
	int m=A->m;

	int N=i2-i1+1;
	int M=j2-j1+1;
	
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			*(A->pointer+m*(i+i1)+j+j1)=*(B->pointer+M*i+j);
		}
	}
}