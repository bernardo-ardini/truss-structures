void _lu_serial(struct matrix *A);
void _lu(struct matrix *A);
void forward_solve_matrix(struct matrix *L,struct matrix *B,struct matrix *X);
void backward_solve_matrix(struct matrix *U,struct matrix *B,struct matrix *X);
void _linsolve(struct matrix *A,struct matrix *b);

void _lu_serial(struct matrix *A)
{
	int n=A->m;
	
	for(int k=0;k<n-1;k++)
	{
		for(int i=k+1;i<n;i++)
		{
			*(A->pointer+i*n+k)=*(A->pointer+i*n+k)/(*(A->pointer+k*n+k));
		}
		
		for(int j=k+1;j<n;j++)
		{
			for(int i=k+1;i<n;i++)
			{
				*(A->pointer+i*n+j)=*(A->pointer+i*n+j)-*(A->pointer+i*n+k)**(A->pointer+k*n+j);
			}
		}
	}
}

void _lu(struct matrix *A)
{
	int N=70;
	int n=A->m;
	
	if(n<=N)
	{
		_lu_serial(A);
		return;
	}
	
	struct matrix A00,A01,A10,A11;
	struct matrix L10;
	struct matrix U01;
	struct matrix B;
	
	readsubm(A,&A00,0,N-1,0,N-1);
	_lu_serial(&A00);
	
	readsubm(A,&A01,0,N-1,N,n-1);
	forward_solve_matrix(&A00,&A01,&U01);
	clearm(&A01);
	
	readsubm(A,&A10,N,n-1,0,N-1);
	backward_solve_matrix(&A00,&A10,&L10);
	clearm(&A10);
	
	readsubm(A,&A11,N,n-1,N,n-1);
	matmul(&L10,&U01,&B);
	_difference(&A11,&B);
	
	_lu(&A11);
	
	_writesubm(A,&A00,0,N-1,0,N-1);
	_writesubm(A,&U01,0,N-1,N,n-1);
	_writesubm(A,&L10,N,n-1,0,N-1);
	_writesubm(A,&A11,N,n-1,N,n-1);
	
	clearm(&A00);
	clearm(&U01);
	clearm(&L10);
	clearm(&A11);
}

double readL(struct matrix *A,int i,int j)
{
	if(i==j) return 1;
	else return *(A->pointer+i*(A->m)+j);
}

void forward_solve_matrix(struct matrix *L,struct matrix *B,struct matrix *X)
{
	int n=L->n;
	int M=B->m;
	
	loadm(X,n,M);
	
	#pragma omp parallel for
	for(int k=0;k<M;k++)
	{
		for(int i=0;i<n;i++)
		{
			*(X->pointer+i*M+k)=*(B->pointer+i*M+k);
			for(int j=0;j<=i-1;j++)
			{
				*(X->pointer+i*M+k)-=*(X->pointer+j*M+k)**(L->pointer+i*n+j);
			}
		}
	}
}

void backward_solve_matrix(struct matrix *U,struct matrix *B,struct matrix *X)
{
	int m=U->m;
	int N=B->n;
	
	loadm(X,N,m);
	
	#pragma omp parallel for
	for(int k=0;k<N;k++)
	{
		for(int i=0;i<m;i++)
		{
			*(X->pointer+i+k*m)=*(B->pointer+i+k*m)/(*(U->pointer+i*m+i));
			for(int j=0;j<=i-1;j++)
			{
				*(X->pointer+i+k*m)-=*(X->pointer+j+k*m)**(U->pointer+j*m+i)/(*(U->pointer+i*m+i));
			}
		}
	}
}

void _forward_solve(struct matrix *L,struct matrix *b)
{
	int n=L->n;
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<=i-1;j++)
		{
			*(b->pointer+i)-=*(b->pointer+j)**(L->pointer+i*n+j);
		}
	}
}

void _backward_solve(struct matrix *U,struct matrix *b)
{
	int n=U->n;
	
	for(int i=n-1;i>=0;i--)
	{
		*(b->pointer+i)/=*(U->pointer+i*n+i);
		for(int j=i+1;j<n;j++)
		{
			*(b->pointer+i)-=*(b->pointer+j)**(U->pointer+i*n+j)/(*(U->pointer+i*n+i));
		}
	}
}

void _linsolve(struct matrix *A,struct matrix *b)
{
	_lu_serial(A);
	_forward_solve(A,b);
	_backward_solve(A,b);
}