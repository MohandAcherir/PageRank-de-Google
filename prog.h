typedef struct Liste {
  int sommet; // Numéro du sommet
  double p; // Proba associée
  struct Liste* suiv; // Suivant dans la liste chainée
} Liste;

void produitVecteurMatrice(Liste **listeSommet, double *vecteur, int taille);
double norme(double *vecteur, int taille);
void printListe(Liste **listeColonnes, int taille);
void add_sommet(Liste **listeColonnes, int NumLigne, int parsedSommet, double parsedProba);
void parse(int *NbSommets, Liste ***listeColonnes, double **F, const char *nomFichier);
void freeListe(Liste *colonne);
void freelisteColonnes(Liste **listeColonnes, int taille);
double prod_scalaire(double *x, double *y, int taille);
double prod_scalaire2(double *x, double alpha, int taille);
void calculxG(Liste **P, double *x, double *f, double alpha, int taille);
double absol(double x);
void normaliser(double *v, int taille);
void copyVecteur(double *dst, double *src, int taille);
int* reorganise (double* x, int NbSommet);
Liste* getHead(Liste* n);
Liste* copySimpleListe(Liste *src);
Liste** copyListe(Liste **src, int taille);
Liste** addSommets(Liste **M, double* x, int NbSommet, double croissance, int degre, double** F);
void converger(Liste **listeColonnes, double *vecteur, double *F, int taille, double *alphaArr, double *precision);
void convergerOnce(Liste **listeColonnes, double *vecteur, double *F, int taille, double alpha, double precision, FILE *fp);


