#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "prog.h"

/**
 * Effectue le produit scalaire entre x et y.
*/
double prod_scalaire(double *x, double *y, int taille) {
  double sum = 0.0;
  for (int i = 0; i < taille; i++) {
    sum += x[i] * y[i];
  }
  return sum;
}

/**
 * Effectue le produit x * (1/N) * (1-α) * e^t * e
*/
double prod_scalaire2(double *x, double alpha, int taille) {
  double sum = 0.0;
  for (int i = 0; i < taille; i++) {
    sum += x[i];
  }
  return (sum * (1-alpha)) / (double) taille;
}

/**
 * Met à jour le vecteur x en le multipliant par la matrice G.
*/
void calculxG(Liste **P, double *x, double *f, double alpha, int taille) {
  double prod1 = prod_scalaire(x, f, taille) * (alpha / (double) taille); // α(1/N)(x f^t)
  produitVecteurMatrice(P, x, taille);
  for (int i = 0; i < taille; i++) {
    x[i] *= alpha;
    // Après cette étape x contient (alpha * x * P), on ajoute donc (1 − α)(1/N) + α(1/N)(x f^t)
    x[i] += (1 - alpha) / taille;
    x[i] += prod1;
  }
}


double absol(double x) {
  if (x < 0) {
    return -x;
  }
  return x;
}

/**
 * Renvoie la norme 1 d'un vecteur.
*/
double norme(double *vecteur, int taille) {
  double sum = 0.0;
  for (int i = 0; i < taille; i++) {
    sum += absol(vecteur[i]);
  }
  return sum;
}

/**
 * Print toutes les colonnes de la matrice.
*/
void printListe(Liste **listeColonnes, int taille) {
  Liste *tete;
  for (int i = 0; i < taille; i++) {
    printf("Colonne %d\n", i + 1);
    tete = listeColonnes[i];
    if (tete != NULL) {
      printf("%d %f\n", tete->sommet + 1, tete->p);
      tete = tete->suiv;
      while (tete != NULL) {
        printf("%d %f\n", tete->sommet + 1, tete->p);
        tete = tete->suiv;
      }
    }
  }
}

/**
 * Ajoute un sommet à une des colonnes de la matrice
*/
void add_sommet(Liste **listeColonnes, int NumLigne, int parsedSommet, double parsedProba) {
  Liste *new = malloc(sizeof(Liste));
  new->sommet = NumLigne;
  new->p = parsedProba;
  if (listeColonnes[parsedSommet] == NULL) {
    listeColonnes[parsedSommet] = new;
  } else {
    new->suiv = listeColonnes[parsedSommet];
    listeColonnes[parsedSommet] = new;
  }
}

/**
 * Parse le fichier, initialise en mémoire et en valeur la matrice et le vecteur F.
*/
void parse(int *NbSommets, Liste ***listeColonnes, double **F, const char *nomFichier) {
  FILE *fd = fopen(nomFichier, "r");
  if (fd == NULL) {
    printf("Impossible d'ouvrir le fichier %s", nomFichier);
    exit(1);
  }
  int NbArcs;
  fscanf(fd, "%d\n%d\n", NbSommets, &NbArcs);
  printf("NbSommets: %d NbArcs: %d\n", *NbSommets, NbArcs);
  *listeColonnes = malloc(*NbSommets * sizeof(Liste));
  *F = malloc(*NbSommets * sizeof(double));
  printf("Malloc reussi\n");
  int NumLigne, NbArcsLigne = 0;
  int parsedSommet;
  double parsedProba;

  for (int i = 0; i < *NbSommets; i++) {
    fscanf(fd, "%d %d", &NumLigne, &NbArcsLigne);
    (*F)[i] = (NbArcsLigne == 0) ? 1.0 : 0.0; // Si un sommet est dans le fichier, sa proba est > 0
    for (int k = 0; k < NbArcsLigne; k++) {
      fscanf(fd, " %d %lf", &parsedSommet, &parsedProba);
      add_sommet(*listeColonnes, NumLigne - 1, parsedSommet - 1, parsedProba);
    }
    fgetc(fd); // Consomme le retour à la ligne
  }
  printf("Parse reussi\n");
  fclose(fd);
}

/**
 * Met à jour le vecteur v en le multipliant par la matrice.
*/
void produitVecteurMatrice(Liste **listeSommet, double *vecteur, int taille) {
  double *newVector = malloc(sizeof(double) * taille);
  double sum;
  Liste *curListe;
  for (int i = 0; i < taille; i++) {
    sum = 0.;
    curListe = listeSommet[i];
    while (curListe != NULL) {
      sum += curListe->p * vecteur[curListe->sommet];
      curListe = curListe->suiv;
    }
    newVector[i] = sum;
  }
  for (int i = 0; i < taille; i++) {
    vecteur[i] = newVector[i];
  }
  free(newVector);  
}

/**
 * Normalise le vecteur (divise par sa norme)
*/
void normaliser(double *v, int taille) {
  double somme = 0.;
  for (int i = 0; i < taille; i++) {
    somme += absol(v[i]);
  }
  for (int i = 0; i < taille; i++) {
    v[i] /= somme;
  }
}


/**
 * Libere une colonne.
*/
void freeListe(Liste *colonne) {
  Liste *tmp;
  while (colonne != NULL) {
    tmp = colonne;
    colonne = colonne->suiv;
    free(tmp);
  }
}

/**
 * Libere toutes les colonnes de la matrice.
*/
void freelisteColonnes(Liste **listeColonnes, int taille) {
  for (int i = 0; i < taille; i++) {
    freeListe(listeColonnes[i]);
  }
  free(listeColonnes);
}

/**
 * Copie le vecteur src vers le vecteur dst
*/
void copyVecteur(double *dst, double *src, int taille) {
  for (int i = 0; i < taille; i++) {
    dst[i] = src[i];
  }
}

double distance(double *vecteur1, double *vecteur2, int taille) {
  double sum = 0.;
  for (int i = 0; i < taille; i++) {
    sum += absol(vecteur1[i] - vecteur2[i]);
  }
  return sum;
}


/**
 * Renvoie les 1% meilleurs sommets en ordre decroissant selon leurs notes
*/
int* reorganise (double* x, int NbSommet) {
	double* X = malloc(sizeof(double)*NbSommet);
	copyVecteur(X, x, NbSommet);
  int nbMeilleurs = (int) (NbSommet * 0.01);

	int* arr = (int*) malloc(sizeof(int) * nbMeilleurs);
	
  // Selectionner les 1% meilleures notes
	for(int i = 0 ; i < nbMeilleurs ; i++) {
		int max = i;
		for(int j = 0 ; j < NbSommet ; j++) {
			if(X[max] < X[j]) {
				max = j;
			}
		}
		X[max] = -1.0;
		arr[i] = max;
	}
	
	free(X);
	return arr;
}


/**
 * Renvoie le dernier element de la liste
*/
Liste* getHead(Liste* n) {
	if(n == NULL) {
		return NULL;
	}
	
	while(n->suiv != NULL) {
		n = n->suiv;
	}
	return n;
}

Liste* copySimpleListe(Liste *src) {
  if (src == NULL) {
    return NULL;
  }
  Liste *new = malloc(sizeof(Liste));
  new->p = src->p;
  new->sommet = src->sommet;
  new->suiv = copySimpleListe(src->suiv);
  return new;
}

Liste** copyListe(Liste **src, int taille) {
  Liste **dst = malloc(taille * sizeof(Liste));
  for (int i = 0; i < taille; i++) {
    dst[i] = copySimpleListe(src[i]);
  }
  return dst;
}

/**
 * Ajoute des sommets et arcs au graphe initial
*/
Liste** addSommets(Liste **M, double* x, int NbSommet, double croissance, int degre, double** F)
{
  Liste **newM = copyListe(M, NbSommet);
  int nbMeilleurs = (int) (NbSommet * 0.01);


	int n = NbSommet*(1+croissance);
	double* newF = malloc(sizeof(double) * n);
	
	for(int i = 0 ; i < n ; i++) {
		if(i < NbSommet) {
			newF[i] = (*F)[i];
		}
		else {
			newF[i] = 1.0;
		}
	}
	srand(time(NULL));
	int extraSommet = NbSommet*croissance; // croissance : Pourcentage de croissance
	int* select = reorganise(x, NbSommet);

	// Pour chaque nouveau sommet s
	int countArcs = 0, countRands = 0;
  int randIndice;
  int *randIndiceTire = malloc(sizeof(int) * degre);
  int currentIndice;
	for (int i = 0; (i < extraSommet); i++)
	{
		int s = NbSommet + i; // Le numero de s
		countRands++;
		//selectionner les nRand meilleurs sommets
		for (int j = 0; (j < degre); j++)
		{
      randIndice = rand() % nbMeilleurs;
      currentIndice = j - 1;
      while (currentIndice >= 0) {
        if (randIndiceTire[currentIndice] == randIndice) {
          currentIndice = j-1;
          randIndice = rand() % nbMeilleurs;
        } else {
          currentIndice--;
        }
      }
      randIndiceTire[j] = randIndice;
      //printf("%d\n", randIndice);
      Liste* new = (Liste*) malloc(sizeof(Liste));
			new->sommet = s;
			new->p = 1.0 / (double) degre;
			new->suiv = NULL;
      // TODO : mémoriser les valurs déjà tirées
			Liste* head = getHead(newM[select[randIndice]]);
			if(head == NULL) {
				newM[select[randIndice]] = new;
			}
			else {
				head->suiv = new;
			}
			countArcs++;
			newF[s] = 0.0;
		}
	}
	
	free(select);	
	printf("Arcs utilises %d, Randoms: %d\n", countArcs, countRands);
	// Nouvelle liste
	Liste** M1 = malloc(sizeof(Liste*) * n);
	for(int i = 0 ; i < n ; i++) {
		if(i < NbSommet) {
			M1[i] = newM[i];
		}
		else {
			M1[i] = NULL;
		}
	}
	free(*F);
	(*F) = newF;
	return M1;
}

#define NB_ALPHA 6
#define NB_PRECISION 7

void converger(Liste **listeColonnes, double *vecteur, double *F, int taille, double *alphaArr, double *precision) {
  int i = 0;
  double *oldVecteur = malloc(taille * sizeof(double));
  for (int k = 0; k < NB_ALPHA-5; k++) {
    for (int h = 0; h < taille; h++) {
      vecteur[h] = 1.0 / (double) taille;
    }
    for (int j = 0; j < NB_PRECISION; j++) {
      do {
        copyVecteur(oldVecteur, vecteur, taille);
        calculxG(listeColonnes, vecteur, F, alphaArr[k], taille);
        normaliser(vecteur, taille);
        i++;
      } while (distance(oldVecteur, vecteur, taille) > precision[j]);
      printf("Convergence en %d iterations pour alpha = %lf précision :10e-%d\n", i, alphaArr[k], j + 3);
    }
    i = 0;
    printf("\n\n");
  }
  free(oldVecteur);
}


void convergerOnce(Liste **listeColonnes, double *vecteur, double *F, int taille, double alpha, double precision, FILE *fp) {
  double *oldVecteur = malloc(taille * sizeof(double));
  int i = 0;
  do {
    copyVecteur(oldVecteur, vecteur, taille);
    calculxG(listeColonnes, vecteur, F, alpha, taille);
    normaliser(vecteur, taille);
    i++;
  } while (distance(oldVecteur, vecteur, taille) > precision);
  fprintf(fp, "%d itérations\n", i);
  fflush(fp);
}


int main(int argc, char const *argv[])
{
  /**
   * Stratégie :
   * - Faire converger avec e/n une premiere fois
   * - Récupérer le vecteur, ajouter des sommets et arcs
   * - Sur la nouvelle matrice, on fait tourner deux algo en parallèle :
   *   - Iteration avec e/n
   *   - Iteration avec note = 0 pour les nouveaux sommets
   * - Comparer les résultats obtenus
  */
  if(argc != 5) {
  	printf("Usage : ./prog croissance nbArcs <nomFichierEntrée> <nomFichierSortie>\n");
  	exit(-1);
  }
  

  FILE *fp = fopen(argv[4], "w+");
  fprintf(fp, "%s\ncroissance : %s\ndegré nouveaux sommets : %s\n\n\n", argv[3], argv[1], argv[2]);
  if (fp == NULL) {
    printf("Impossible d'ouvrir le fichier %s", argv[4]);
  }



  double alphaArr[NB_ALPHA] = {0.5,0.7,0.85,0.9,0.99,0.999};
  double precision = 1e-9;

  int taille;
  Liste **M1;
  double *F;

  parse(&taille, &M1, &F, argv[3]);


  double *F1 = malloc(sizeof(double) * taille);
  double *F2 = malloc(sizeof(double) * taille);
  for (int i = 0; i < taille; i++) {
    F1[i] = F[i];
  }

  
    
  double croissance = strtod(argv[1], NULL);  //Croissance
  int degre = atoi(argv[2]);  //Arcs supplementaires
  int newTaille = taille * (1 + croissance);

  double *vecteur = malloc(taille * sizeof(double));
  // Méthode initiale avec e/n
  double *vecteurBasique = malloc(sizeof(double) * newTaille);
  // Méthode attachement préferentiel
  double *vecteurPref = malloc(sizeof(double) * newTaille);

  Liste** M2;
  for (int a = 0; a < NB_ALPHA; a++) {  
    fprintf(fp, "Alpha = %lf\n", alphaArr[a]);
    for (int i = 0; i < taille; i++) {
      vecteur[i] = 1.0 / (double) taille;
    }
    fprintf(fp, "Calcul de note :");
    convergerOnce(M1, vecteur, F1, taille, alphaArr[a], precision, fp);
    copyVecteur(F2, F1, taille);
    M2 = addSommets(M1, vecteur, taille, croissance, degre, &F2);
    for (int i = 0; i < newTaille; i++) {
      vecteurBasique[i] = 1.0 / (double) newTaille;
      if (i < taille) {
        vecteurPref[i] = vecteur[i];
      } else {
        vecteurPref[i] = 0.0;
      }
    }
    fprintf(fp, "basique :");
    convergerOnce(M2, vecteurBasique, F2, newTaille, alphaArr[a], precision, fp);
    fprintf(fp, "pref :");
    convergerOnce(M2, vecteurPref, F2, newTaille, alphaArr[a], precision, fp);
    fprintf(fp, "\n\n");
    freelisteColonnes(M2, newTaille);
  }
  
  freelisteColonnes(M1, taille);
  free(F);
  free(F1);
  free(F2);
  free(vecteur);
  free(vecteurPref);
  free(vecteurBasique);
  return 0;  	
}
