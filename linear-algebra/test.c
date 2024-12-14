#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <assert.h>
#include"matrix.c"
#include"lu.c"

int main()
{	
	struct matrix A;

	loadm(&A,5,5);
	setm(&A,0,0,17);
	setm(&A,0,1,24);
	setm(&A,0,2,1);
	setm(&A,0,3,8);
	setm(&A,0,4,15);
	setm(&A,1,0,23);
	setm(&A,1,1,5);
	setm(&A,1,2,7);
	setm(&A,1,3,14);
	setm(&A,1,4,16);
	setm(&A,2,0,4);
	setm(&A,2,1,6);
	setm(&A,2,2,13);
	setm(&A,2,3,20);
	setm(&A,2,4,22);
	setm(&A,3,0,10);
	setm(&A,3,1,12);
	setm(&A,3,2,19);
	setm(&A,3,3,21);
	setm(&A,3,4,3);
	setm(&A,4,0,11);
	setm(&A,4,1,18);
	setm(&A,4,2,25);
	setm(&A,4,3,2);
	setm(&A,4,4,9);
	
	// test linsolve
	
	struct matrix b,x;
	
	ones(&x,5,1);
	matmul(&A,&x,&b);
	_linsolve(&A,&b);
	printm(&b);
	
	// test LU parallel

	/*

    printm(&A);
    printf("\n");
    _lu(&A);
    printm(&A);
    
    */
	
	// test LU serial
	
	/*
    printm(&A);
    printf("\n");
    _lu_serial(&A);
    printm(&A);
    */
    
    // test readsubm
    
    /*
    
    struct matrix B;
    readsubm(&A,&B,3,4,1,2);
    printm(&A);
    printf("\n");
    printm(&B);
    
    */
    
    // test forward solve matrix
    
    /*
    loadm(&A,2,2);
    setm(&A,0,0,1);
    setm(&A,0,1,0);
    setm(&A,1,0,2);
    setm(&A,1,1,1);
    printm(&A);
    printf("\n");
    
    struct matrix B;
    loadm(&B,2,3);
    setm(&B,0,0,1);
    setm(&B,1,0,2);
    setm(&B,0,1,3);
    setm(&B,1,1,4);
    setm(&B,0,2,5);
    setm(&B,1,2,6);
    printm(&B);
    printf("\n");
    
    struct matrix X;
    forward_solve_matrix(&A,&B,&X);
    
    printm(&X);
    
    */
    
    // test backward solve matrix
    
    /*
    loadm(&A,2,2);
    setm(&A,0,0,1);
    setm(&A,0,1,2);
    setm(&A,1,0,0);
    setm(&A,1,1,1);
    printm(&A);
    printf("\n");
    
    struct matrix B;
    loadm(&B,2,3);
    setm(&B,0,0,1);
    setm(&B,1,0,2);
    setm(&B,0,1,3);
    setm(&B,1,1,4);
    setm(&B,0,2,5);
    setm(&B,1,2,6);
    printm(&B);
    printf("\n");
    
    struct matrix X;
    backward_solve_matrix(&A,&B,&X);
    
    printm(&X);
    */
}