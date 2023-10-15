#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double a = 0.5;

struct Liste{
	int sommet;
	double p;
	struct Liste* suiv;
};

typedef struct Liste Liste;

Liste* getHead(Liste* n) {
	if(n == NULL) {
		return NULL;
	}
	
	while(n->suiv != NULL) {
		n = n->suiv;
	}
	return n;
}

double dot(int f[], double x[], int n, int flag) {
	double ret = 0.0;
	if(flag == 0) {
		for(int i = 0 ; i < n ; i++) {
			ret += f[i]*x[i];
		}
	}
	if(flag == 1){
		for(int i = 0 ; i < n ; i++) {
			ret += x[i]*((double)(1.0/n));
		}
	}
	return ret;
}

void mult(Liste* n, int i, double* v, double* res) {
	if(n != NULL) {
		(*res) += v[(n->sommet)-1] * (n->p);
		
		if(n->suiv != NULL) {
			mult(n->suiv, i+1, v, res);
		}
	}
}


void calcul_xG(Liste* M[], double v[], double u[], int F[], int NbSommet) {
	// Calcul de xG1
	double* G1 = malloc(sizeof(double)*NbSommet);
	for(int i = 1 ; i < NbSommet+1 ; i++) {
		G1[i-1] = 0.0;
		mult(M[i], i, v, &G1[i-1]);
		G1[i-1] *= a;
	}
	// Calcul de xG2
	double* G2 = malloc(sizeof(double)*NbSommet);
	double b = dot(F, v, NbSommet, 0);
	for(int i = 1 ; i < NbSommet+1 ; i++) {
		G2[i-1] = b;
		G2[i-1] *= a;
	}
	// Calcul de xG3
	double* G3  =  malloc(sizeof(double)*NbSommet);
	b = dot(F, v, NbSommet, 1);
	for(int i = 1 ; i < NbSommet+1 ; i++) {
		G3[i-1] = b;
		G3[i-1] *= (1-a);
	}
	// Calcul de xG
	for(int i = 1 ; i < NbSommet+1 ; i++) {
		u[i-1] = G1[i-1]+G2[i-1]+G3[i-1];
	}
	
	free(G1);
	free(G2);
	free(G3);
}

double norme(double* u, double* v, int n) {
	double dist = 0.0;
	
	for(int i = 0 ; i < n ; i++) {
		dist+= (u[i]-v[i])*(u[i]-v[i]);
	}
	return sqrt(dist);
}

void affiche(Liste* n) {
	if(n != NULL) {
		printf("(%d, %.1f) -> ", n->sommet, n->p);
		
		if(n->suiv != NULL) {
			affiche(n->suiv);
		}
		else{
			printf("NULL\n");
		}
	}
}


void freeListe(Liste* n) {
	Liste* s = NULL;
	if(n != NULL) {
		s = n->suiv;
		free(n);
	}	
	
	if(s != NULL) {
		freeListe(s);
	}	
}

void copy(double* u, double* v, int n) {
	for(int i = 0 ; i < n ; i++) {
		v[i] = u[i];
		u[i] = 0.0;
	}
}


// calculer le vecteur ligne f
void calcul_f(Liste* M[], int* tab_f, double* v, double* u, double* zeros, int LENGTH) {
	for(int i = 1 ; i < LENGTH+1 ; i++) {
		double somme = 0;
		tab_f[i-1] = 0;
		v[i-1] = (double)1.0 / LENGTH;
		u[i-1] = 0.0;
		zeros[i-1] = 0.0;
		/*
		for(int j = 1 ; j < LENGTH+1 ; j++) {
			//Liste* curr = search(M[j], i);
			Liste *curr = M[j];
			
			while( (curr != NULL) && (curr->sommet != i) )
			{	curr = curr->suiv;	}
			
			if(curr != NULL) {
				somme += curr->p;
			}
		}
		if(somme == 0.0) {
			tab_f[i-1] = 1;
		}*/
	}
}

void normaliser(double* v, double n, int NbSommet) {
	for(int i = 0 ; i < NbSommet ; i++) {
		v[i] /= n;
	}
}

int main(int argc, char** argv) {
	if(argc != 3) {
		printf("Usage : ./td file.txt nb\n");
		exit(-1);
	}
	FILE* f = fopen(argv[1], "r");
	int NbSommet = 0, NbArcs = 0;
	
	fscanf (f, "%d", &NbSommet);
	fscanf (f, "%d", &NbArcs);
	
	printf("Nb Sommets :%d , NbArcs: %d\n", NbSommet, NbArcs);
	Liste** M = malloc(sizeof(Liste*) * (NbSommet + 1));
	
	double i = 0.0;
	int t = 0;
	double tarc = 0;
	double succ;
	int pred;
	double success = 0.0;
	int* F = malloc(sizeof(int)*NbSommet);
	
	for(int k = 0 ; k < NbSommet+1 ; k++) {
		M[k] =  NULL;
		
		if(k < NbSommet)
			F[k] = 1;
	}
	
	while(!feof(f)) {
		//printf("\nt : %d", t);
		if(t == 0) {
			fscanf (f, "%lf", &i);
		}
		else if(t == 1) {
			fscanf (f, "%lf", &tarc);
		}
		else if(t > 1) {
			if(t%2 == 0) {
				fscanf(f, "%lf", &succ);
				pred = (int)succ;
				
				if(getHead(M[pred]) == NULL) {
					M[pred] = (Liste*)malloc(sizeof(Liste));
					M[pred]->suiv = NULL;
				}
				else {
					getHead(M[pred])->suiv = (Liste*)malloc(sizeof(Liste));
					getHead(M[pred])->suiv = NULL;
				}
				F[((int)pred)-1] = 0;
			}
			if(t%2 == 1) {
				fscanf (f, "%lf", &success);
				
				getHead(M[pred])->sommet = i;
				getHead(M[pred])->p = success;
			}
		}
		
		t += 1;
		if(t >= (2 + 2*tarc)) {
			t = 0;
		}
	}
	fclose(f);
	
	double* v = malloc(sizeof(double)*NbSommet);	
	double* u = malloc(sizeof(double)*NbSommet);
	double* zeros = malloc(sizeof(double)*NbSommet);
	
	calcul_f(M, F, v, u, zeros, NbSommet);
	
	printf("Reading Finished!\n");
	
	int counter = 0;
	int ep = 4;
	do {
		if(counter > 0) {
			//normaliser(u, norme(u, zeros, NbSommet), NbSommet);
			copy(u, v, NbSommet);
		}	

		calcul_xG(M, v, u, F, NbSommet);
		
		counter += 1;
		//printf("Counter : %d\n", counter);
		if(norme(u, v, NbSommet) < ((double)1.0/pow(10, ep))) {
			printf("norme < 10^(-%d) after %d iterations", ep, counter);
			ep++;
		}
	}while(norme(u, v, NbSommet) > ((double)1.0/pow(10, atoi(argv[2]))));
	
	printf("Converges after looping ....%d times\n", counter	);
	
	/*for(int i = 0 ; i < NbSommet + 1 ; i++) {
		if(i == 0) printf("\n(");
		if(i == NbSommet-1) {
			printf("%0.2f)\n", v[i]);
		}else {
			if(i < NbSommet) {
				printf("%0.2f, ", v[i]);
			}	
		}
		//freeListe(M[i]);
	}*/
		
	for(int i = 0 ; i < NbSommet + 1 ; i++) {
		freeListe(M[i]);
	}
	free(M);
	free(v);
	free(u);
	free(F);
	free(zeros);
	
	return 0;
}


//Avec epsilon = 10^(-6)
//graphe 101 : 13 iterations
//graphe 1000: 9  iterations
//graphe 100001: 10 iteraions , prend ~6 secs de temps et 16763 pour 10^(-8), et 13 pour 10^(-7)
//graphe wb-cs-stanford : 14 iterations ~8 secs
//graph Stanford : 14 iterations pour 10^(-3), 29 pour -4, 148 pour -5, 258 pour -6

