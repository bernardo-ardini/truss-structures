void nr(struct matrix *X,struct matrix *x,void (*K)(struct matrix *,struct matrix *),void (*T)(struct matrix *,struct matrix *),struct matrix *F,double tol)
{
	int n=X->n;
	copy(X,x);
	
	struct matrix k,t;
	struct matrix f,r;
	
	(*T)(x,&t);
	(*K)(x,&k);
	
	for(int q=0;q<(F->m);q++)
	{
		if(q>0)
		{
			clearm(&r);
			clearm(&f);
		}
		column(F,q,&f);
		difference(&t,&f,&r);
		
		printf("%2.0f/%2.0f        ",(double)q+1,(double)F->m);
		
		int c=0;
		
		while(norm(&r)>tol*norm(&f))
		{
			c++;
			printf("\b\b\b\b\b\b%2.0f    ",(double)c);
			_product(-1,&k);
			_linsolve(&k,&r);
			_sum(x,&r);
			clearm(&t);
			clearm(&k);
			(*K)(x,&k);
			(*T)(x,&t);
			clearm(&r);
			difference(&t,&f,&r);
		}
		
		printf("\n");
	}
	
	clearm(&t);
	clearm(&k);
	clearm(&f);
	clearm(&r);
}

void nrs(struct matrix *X,struct matrix *x,struct matrix *xs,void (*K)(struct matrix *,struct matrix *),void (*T)(struct matrix *,struct matrix *),struct matrix *F,double tol)
{
	int n=X->n;
	copy(X,x);
	
	struct matrix k,t;
	struct matrix f,r;
	
	(*T)(x,&t);
	(*K)(x,&k);
	
	for(int q=0;q<(F->m);q++)
	{
		if(q>0)
		{
			clearm(&r);
			clearm(&f);
		}
		column(F,q,&f);
		difference(&t,&f,&r);
		
		printf("%2.0f/%2.0f        ",(double)q+1,(double)F->m);
		
		int c=0;
		
		while(norm(&r)>tol*norm(&f))
		{
			c++;
			printf("\b\b\b\b\b\b%2.0f    ",(double)c);
			_product(-1,&k);
			_linsolve(&k,&r);
			_sum(x,&r);
			clearm(&t);
			clearm(&k);
			(*K)(x,&k);
			(*T)(x,&t);
			clearm(&r);
			difference(&t,&f,&r);
		}
		
		copy(x,xs+q);
		
		printf("\n");
	}
	
	clearm(&t);
	clearm(&k);
	clearm(&f);
	clearm(&r);
}