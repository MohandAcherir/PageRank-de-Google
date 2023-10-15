#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double a = 0.7;

// La liste des Alpha
double A[6] = {0.5, 0.7, 0.85, 0.9, 0.99, 0.999};

struct Liste
{
	int sommet;
	double p;
	struct Liste *suiv;
};

typedef struct Liste Liste;

Liste *getHead(Liste *n)
{
	if (n == NULL)
	{
		return NULL;
	}

	while (n->suiv != NULL)
	{
		n = n->suiv;
	}
	return n;
}

/**
 * Renvoie le produit scalaire entre le vecteur v1 et v2.
*/
double produit_scalaire(int v1[], double v2[], int n)
{
	double ret = 0.;
	for (int i = 0; i < n; i++)
	{
		ret += v1[i] * v2[i];
	}
	return ret;
}

/**
 * Renvoie le produit scalaire entre le vecteur v et le vecteur e / n = (1/n,.., 1/n).
*/
double produit_scalaire_e(double v[], int n)
{
	double ret = 0.;
	for (int i = 0; i < n; i++)
	{
		ret += v[i];
	}
	return ret / (double) n;
}

// Une multiplication pour chaque valeur du vecteur d'arrivée (?)
void mult(Liste *n, int i, double *v, double *res)
{
	if (n != NULL)
	{
		(*res) += v[(n->sommet) - 1] * (n->p);

		if (n->suiv != NULL)
		{
			mult(n->suiv, i + 1, v, res);
		}
	}
}

void calcul_xG(Liste *M[], double v[], double u[], int F[], int NbSommet)
{
	// Calcul de xG1
	double *G1 = malloc(sizeof(double) * NbSommet);
	for (int i = 0; i < NbSommet; i++)
	{
		G1[i] = 0.0;
		mult(M[i], i, v, &G1[i]);
		G1[i] *= a;
	}
	// Calcul de xG2
	double *G2 = malloc(sizeof(double) * NbSommet);
	double b = produit_scalaire(F, v, NbSommet);
	for (int i = 0; i < NbSommet; i++)
	{
		G2[i] = b * a;
	}
	// Calcul de xG3
	double *G3 = malloc(sizeof(double) * NbSommet);
	b = produit_scalaire_e(v, NbSommet);
	for (int i = 0; i < NbSommet; i++)
	{
		G3[i] = b * (1 - a);
	}
	// Calcul de xG
	for (int i = 0; i < NbSommet; i++)
	{
		u[i] = G1[i] + G2[i] + G3[i];
	}

	free(G1);
	free(G2);
	free(G3);
}

double val_abs(double x) {
	if (x < 0.) {
		return -x;
	}
	return x;
}

// Distance entre le vecteur u et v selon la norme 1.
double distance(double *u, double *v, int n)
{
	double dist = 0.0;
	for (int i = 0; i < n; i++)
	{
		dist += val_abs(v[i] - u[i]);
	}
	return dist;
}

// Afficher la matrice par colonnes

/*
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
}*/

void freeListe(Liste *n)
{
	Liste *s = NULL;
	if (n != NULL)
	{
		s = n->suiv;
		free(n);
	}

	if (s != NULL)
	{
		freeListe(s);
	}
}

/**
 * Copie le vecteur u dans le vecteur v, le vecteur u est remis à zéro.
*/
void copy_set_zero(double *u, double *v, int n)
{
	for (int i = 0; i < n; i++)
	{
		v[i] = u[i];
		u[i] = 0.0;
	}
}

/**
 * Initialise un vecteur de longueur n avec une valeur donnée.
*/
void init(double *v, double val, int n) {
	for (int i = 0; i < n; i++)
	{
		v[i] = val;
	}
} 

/**
 * Divise le vecteur v par un scalaire x.
*/
void normaliser(double *v, double x, int NbSommet)
{
	for (int i = 0; i < NbSommet; i++)
	{
		v[i] /= x;
	}
}

void parse(FILE *f, int NbSommet, int NbArcs, Liste **M, int *F)
{
	double i = 0.0;
	int t = 0;
	double tarc = 0;
	double succ;
	int pred;
	double success = 0.0;

	for (int k = 0; k < NbSommet + 1; k++)
	{
		M[k] = NULL;

		if (k < NbSommet)
			F[k] = 1;
	}

	while (!feof(f))
	{
		if (t == 0)
		{
			fscanf(f, "%lf", &i);
		}
		else if (t == 1)
		{
			fscanf(f, "%lf", &tarc);
		}
		else if (t > 1)
		{
			if (t % 2 == 0)
			{
				fscanf(f, "%lf", &succ); // parser directement un int %d ?
				pred = (int)succ;

				if (getHead(M[pred]) == NULL)
				{
					M[pred] = (Liste *)malloc(sizeof(Liste));
					M[pred]->suiv = NULL;
				}
				else
				{
					getHead(M[pred])->suiv = (Liste *)malloc(sizeof(Liste));
					getHead(M[pred])->suiv = NULL;
				}
				F[((int)pred) - 1] = 0; // pred déjà un int 
			}
			if (t % 2 == 1)
			{
				fscanf(f, "%lf", &success);

				getHead(M[pred])->sommet = i;
				getHead(M[pred])->p = success;
			}
		}

		t += 1;
		if (t >= (2 + 2 * tarc))
		{
			t = 0;
		}
	}
	fclose(f);
}

void converger(Liste **M, double *u, double *v, double *zeros, int *F, int NbSommet, int nb)
{
	for (int y = 0; y < 6; y++)
	{
		a = A[y];
		printf("\na = %lf\n", a);
		int counter = 0;
		int ep = 3;
		do
		{
			// déplacer ces instructions à la fin de la boucle ?
			if (counter > 0)
			{
				normaliser(u, distance(u, zeros, NbSommet), NbSommet);
				copy_set_zero(u, v, NbSommet);
			}

			// u <- v.G
			calcul_xG(M, v, u, F, NbSommet);

			counter += 1;
			// printf("Counter : %d\n", counter);
			if (distance(u, v, NbSommet) < ((double)1.0 / pow(10, ep)))
			{
				printf("distance < 10^(-%d) after %d iterations\n", ep, counter);
				ep++;
			}
		} while (distance(u, v, NbSommet) > ((double)1.0 / pow(10, nb)));
	}
}

// Ajouter des sommets des Arcs
int addSommets(Liste **M, int NbSommet)
{
	int extraSommet = 4; // Rand between 2% and 10%
	// 1 - Selectionner les ${extraSommet} sommets qui ont la plus grande note
	int sommets[4] = {1, 9, 14, 19};
	size_t n = sizeof(sommets) / sizeof(sommets[0]);
	// 2 - Rajouter a M selon l'etape 1 les nouveau sommets et leurs arcs
	// ici les Arcs sont distribués uniformément sur les sommets séléctionnés.
	for (int i = 1; i <= extraSommet; i++)
	{
		int s = NbSommet + i;
		for (int j = 0; j < n; j++)
		{
			Liste *new = (Liste *)malloc(sizeof(Liste));
			new->sommet = s;
			new->p = 1.0 / (double)n;
			new->suiv = NULL;

			Liste *head = getHead(M[sommets[j]]);
			if (head == NULL)
			{
				M[sommets[j]] = new;
			}
			else
			{
				head->suiv = new;
			}
		}
	}
	return NbSommet + extraSommet;
}

int main(int argc, char **argv)
{
	// Préliminaires
	if (argc != 3)
	{
		printf("Usage : ./td file.txt nombre\n");
		// nombre : faires les convergences de 10^(-3) à 10^(-nombre)
		exit(-1);
	}
	FILE *f = fopen(argv[1], "r");

	if (f == NULL)
	{
		printf("Error, file does not exist!\n");
		exit(-1);
	}

	// Partie Parsing
	int NbSommet, NbArcs;
	fscanf(f, "%d", &NbSommet);
	fscanf(f, "%d", &NbArcs);
	printf("Nb Sommets :%d , NbArcs: %d\n", NbSommet, NbArcs);

	Liste **M = malloc(sizeof(Liste *) * (NbSommet + 1));
	int *F = malloc(sizeof(int) * NbSommet);

	parse(f, NbSommet, NbArcs, M, F);
	printf("Reading Finished!\n");

	// Partie Convergence
	double *v = malloc(sizeof(double) * NbSommet);
	double *u = malloc(sizeof(double) * NbSommet);
	double *zeros = malloc(sizeof(double) * NbSommet);

	// initialiser les vecteur u, v, zeros créés
	init(u, 0., NbSommet);
	init(zeros, 0., NbSommet);
	init(v, 1. / NbSommet, NbSommet);

	// Converger
	converger(M, u, v, zeros, F, NbSommet, atoi(argv[2]));

	// Ajout de sommets, ici 4 pour tester
	int NewNB = addSommets(M, NbSommet);
	// Reinitialiser les vecteur v, u, zeros et F
	// Converger

	for (int i = 0; i < NbSommet + 1; i++)
	{
		freeListe(M[i]);
	}
	free(M);
	free(v);
	free(u);
	free(F);
	free(zeros);

	return 0;
}
