struct structure
{
	int n;
	int nelements;
	struct element *elements;
};

struct structure str,substr; // structure and current substructure
struct matrix X,subX; // referential coordinates of structure and substructure
struct matrix x; // spatial coordinates of structure
struct matrix *xs; // hystory of spatial coordinates of structure
int *A; // indexes of free nodes
int m; // number of free nodes
double Xtol=1e-2; // tolerance to see if two nodes coincides
struct matrix *F; // external forces
int N=10; // number of steps

double g=98; // acceleration of gravity

void T(struct matrix *y,struct matrix *t);
void compute_x(struct matrix *y,struct matrix *x);
void K(struct matrix *y,struct matrix *k);
int node_yet(int a);
void add_element(struct element *e);
void add_substructure();

void T(struct matrix *y,struct matrix *t)
{
	zeros(t,3*m,1);
	struct matrix x;
	compute_x(y,&x);
	#pragma omp parallel
	{
	#pragma omp single
	{
	for(int a=0;a<m;a++)
	{
		#pragma omp task firstprivate(a)
		{
		for(int l=0;l<str.nelements;l++)
		{
			struct element *e;
			e=str.elements+l;
			int s1=s(e,*(A+a));
			if(s1!=0)
			{
				struct matrix tt;
				Te(e,&x,&tt);
				for(int i=0;i<3;i++)
				{
					*(t->pointer+3*a+i)+=s1**(tt.pointer+i);
				}
				clearm(&tt);
			}
		}
		}
	}
	}
	}
	clearm(&x);
}

void K(struct matrix *y,struct matrix *k)
{
	zeros(k,3*m,3*m);
	struct matrix x;
	compute_x(y,&x);
	#pragma omp parallel
	{
	#pragma omp single
	{
	for(int a=0;a<m;a++)
	{
		for(int b=0;b<m;b++)
		{
			#pragma omp task firstprivate(a,b)
			{
			for(int l=0;l<str.nelements;l++)
			{
				struct element *e;
				e=str.elements+l;
				int s1=s(e,*(A+a));
				int s2=s(e,*(A+b));
				if(s1!=0 && s2!=0)
				{
					struct matrix tt;
					Ke(e,&x,&tt);
					for(int i=0;i<3;i++)
					{
						for(int j=0;j<3;j++)
						{
							*(k->pointer+3*3*m*a+3*b+3*m*i+j)+=s1*s2**(tt.pointer+i*3+j);
						}
					}
					clearm(&tt);
				}
			}
			}
		}
	}
	}
	}
	clearm(&x);
}

void compute_x(struct matrix *y,struct matrix *x)
{
	copy(&X,x);
	
	for(int a=0;a<m;a++)
	{
		for(int i=0;i<3;i++)
		{
			int l=x->m;
			*((x->pointer)+*(A+a)+l*i)=*(y->pointer+3*a+i);
		}
	}
}

void compute_y(struct matrix *x,struct matrix *y)
{
	loadm(y,3*m,1);
	
	for(int a=0;a<m;a++)
	{
		for(int i=0;i<3;i++)
		{
			int l=x->m;
			*(y->pointer+3*a+i)=*((x->pointer)+*(A+a)+l*i);
		}
	}
}

void compute_G(struct matrix *F,int N,struct matrix *G)
{
	loadm(G,3*m,N);
	
	for(int q=0;q<N;q++)
	{
		struct matrix *f;
		f=F+q;
		
		for(int a=0;a<m;a++)
		{
			for(int i=0;i<3;i++)
			{
				int l=f->m;
				*(G->pointer+N*i+3*N*a+q)=*(f->pointer+*(A+a)+l*i);
			}
		}
	}
}

void equilibrium(double tol)
{
	struct matrix y,Y,G;
	struct matrix *ys;
	
	compute_y(&X,&Y);
	ys=(struct matrix *)malloc((N+1)*sizeof(struct matrix));
	
	compute_G(F,N+1,&G);
	
	nrs(&Y,&y,ys,K,T,&G,tol);
	
	compute_x(&y,&x);
	
	xs=(struct matrix *)malloc((N+1)*sizeof(struct matrix));
	for(int q=0;q<=N;q++)
	{
		compute_x(ys+q,xs+q);
	}
	
	clearm(&y);
	clearm(&Y);
	clearm(&G);
	free(ys);
}

int node_yet(int a)
{
	struct matrix x1;
	column(&subX,a,&x1);
	
	for(int b=0;b<str.n;b++)
	{
		struct matrix x2;
		column(&X,b,&x2);
		
		_difference(&x2,&x1);
		
		if(norm(&x2)<Xtol) return b;
		
		clearm(&x2);
	}
	
	clearm(&x1);
	
	return -1;
}

int get_index_node(struct matrix *x)
{
	for(int b=0;b<str.n;b++)
	{
		struct matrix x2;
		column(&X,b,&x2);
		
		_difference(&x2,x);
		
		if(norm(&x2)<Xtol) return b;
		
		clearm(&x2);
	}
	
	return -1;
}

void free_new_node(int a)
{	
	int *B;
	B=(int*)malloc((m+1)*sizeof(int));
	
	for(int b=0;b<m;b++)
	{
		*(B+b)=*(A+b);
	}
	
	*(B+m)=a;
	if(m>0) free(A);
	A=B;
	m++;
}

void fix_node(struct matrix *x)
{
	int a=get_index_node(x);
	
	if(a==-1) return;
	
	int *B;
	B=(int*)malloc((m-1)*sizeof(int));
	
	int found=0;
	
	for(int b=0;b<m;b++)
	{
		if(*(A+b)!=a)
		{
			*(B+b-found)=*(A+b);
		}
		else
		{
			found=1;
		}
	}
	
	if(m>0) free(A);
	A=B;
	if(found==1) m--;
}

void add_element(struct element *e)
{
	struct element *els;
	els=(struct element*)malloc((str.nelements+1)*sizeof(struct element));
	
	for(int l=0;l<str.nelements;l++)
	{
		copye(str.elements+l,els+l);
	}
	copye(e,els+str.nelements);
	
	if(str.nelements>0) free(str.elements);
	
	str.elements=els;
	str.nelements++;
}

void add_substructure()
{
	for(int l=0;l<substr.nelements;l++)
	{
		struct element *e;
		e=substr.elements+l;
		
		int a=node_yet(e->e1);
		int b=node_yet(e->e2);
		
		int element_yet=0;
		
		if(a==-1)
		{
			struct matrix c;
			column(&subX,e->e1,&c);
			_append_column(&X,&c);
			clearm(&c);
			a=str.n;
			free_new_node(a);
			str.n++;
		}
			
		if(b==-1)
		{
			struct matrix c;
			column(&subX,e->e2,&c);
			_append_column(&X,&c);
			clearm(&c);
			b=str.n;
			free_new_node(b);
			str.n++;
		}
		
		e->e1=a;
		e->e2=b;
		
		for(int m=0;m<str.nelements;m++)
		{
			struct element *f;
			f=str.elements+m;
			
			if((e->e1==f->e1 && e->e2==f->e2) || (e->e1==f->e2 && e->e2==f->e1)) element_yet=1;
		}
		
		if(element_yet==0)
		{
			add_element(e);
		}
	}
}

void inizialize_structure()
{
	X.n=3;
	X.m=0;
	str.n=0;
	str.nelements=0;
	m=0;
}

void inizialize_substructure()
{
	subX.n=3;
	subX.m=0;
	substr.n=0;
	substr.nelements=0;
}

void clear_substructure()
{
	if(substr.nelements>0) free(substr.elements);
	if(substr.n>0) clearm(&subX);
	substr.nelements=0;
	substr.n=0;
}

void export_elements(char *filename)
{
	FILE *file;
	file=fopen(filename,"w");
	
	for(int l=0;l<str.nelements;l++)
	{
		struct element *e;
		e=str.elements+l;
		fprintf(file,"%d\t%d\n",e->e1,e->e2);
	}

	fclose(file);
}

void export_story(char *filename)
{
	FILE *file;
	file=fopen(filename,"w");
	fprintf(file,"");
	fclose(file);
	
	for(int q=0;q<=N;q++)
	{
		file=fopen(filename,"a");
		fprintf(file,"# time %d\n",q);
		fclose(file);
		fprintm(filename,xs+q,"a");
		file=fopen(filename,"a");
		fprintf(file,"# end\n");
		fclose(file);
	}
	
	file=fopen(filename,"a");
	fprintf(file,"# END");
	fclose(file);
}

void inizialize_force()
{
	F=(struct matrix *)malloc((N+1)*sizeof(struct matrix));
	
	for(int i=0;i<=N;i++)
	{
		zeros(F+i,3,str.n);
	}
}

void add_weight()
{	
	struct matrix W;
	zeros(&W,3,str.n);
	
	for(int a=0;a<str.n;a++)
	{
		for(int l=0;l<str.nelements;l++)
		{
			struct element *e;
			e=str.elements+l;
			int s1=s(e,*(A+a));
			if(s1!=0)
			{
				*(W.pointer+str.n*2+a)+=-g*e->m/2;
			}
		}
	}
	
	for(int i=0;i<=N;i++)
	{
		for(int a=0;a<str.n;a++)
		{
   			*((F+i)->pointer+2*str.n+a)+=*(W.pointer+str.n*2+a)*(double)i/(double)N;
		}
	}
	
	clearm(&W);
}

void add_referential_force(void (*phi)(struct matrix *,struct matrix *))
{	
	struct matrix W;
	zeros(&W,3,str.n);
	
	for(int a=0;a<str.n;a++)
	{
		struct matrix v,x;
		
		column(&X,a,&x);
		phi(&x,&v);
		
		clearm(&x);
		
		setm(&W,0,a,getm(&v,0,0));
		setm(&W,1,a,getm(&v,1,0));
		setm(&W,2,a,getm(&v,2,0));
		
		clearm(&v);
	}
	
	for(int i=0;i<=N;i++)
	{
		for(int a=0;a<str.n;a++)
		{
			for(int j=0;j<3;j++)
			{
   				*((F+i)->pointer+j*str.n+a)+=*(W.pointer+str.n*j+a)*i/N;
   			}
		}
	}
	
	clearm(&W);
}
