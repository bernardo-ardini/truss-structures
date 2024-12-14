struct element
{
	int e1,e2;
	double A,L;
	double (*tau)(double);
	double (*dtau)(double);
	double m;
};

void Te(struct element *e,struct matrix *x,struct matrix *t);
void Ke(struct element *e,struct matrix *x,struct matrix *k);
int s(struct element *e,int a);

void copye(struct element *e,struct element *f)
{
	f->e1=e->e1;
	f->e2=e->e2;
	f->A=e->A;
	f->L=e->L;
	f->tau=e->tau;
	f->dtau=e->dtau;
	f->m=e->m;
}

void Te(struct element *e,struct matrix *x,struct matrix *t)
{
	struct matrix x1,x2;
	column(x,e->e1,&x1);
	column(x,e->e2,&x2);
	difference(&x2,&x1,t);
	clearm(&x1);
	clearm(&x2);
	double l=norm(t);
	double epsilon=log(l/e->L);
	_product(1/l,t);
	double V=e->A*e->L;
	double tau=(*(e->tau))(epsilon);
	_product(tau*V/l,t);
}

void Ke(struct element *e,struct matrix *x,struct matrix *k)
{
	struct matrix x1,x2,n,J;
	column(x,e->e1,&x1);
	column(x,e->e2,&x2);
	difference(&x2,&x1,&n);
	clearm(&x1);
	clearm(&x2);
	double l=norm(&n);
	double epsilon=log(l/e->L);
	_product(1/l,&n);
	double V=e->A*e->L;
	double tau=(*(e->tau))(epsilon);
	double dtau=(*(e->dtau))(epsilon);
	dyadic(&n,k);
	clearm(&n);
	_product(V/pow(l,2)*(dtau-2*tau),k);
	eye(&J,3);
	_product(tau*V/pow(l,2),&J);
	_sum(k,&J);
	clearm(&J);
}

int s(struct element *e,int a)
{
	if(a==e->e1) return -1;
	else if(a==e->e2) return 1;
	else return 0;
}