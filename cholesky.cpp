/**
 * @file RandriamanantenaManitraLuc.cpp
 * @author your name (randiluc@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-04
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
using namespace std;
void echange(double *a, double *b);
void affichetab(vector<vector<double>> tab, vector<double> b);
void solution(int taille, vector<vector<double>> &tab, vector<double> &b);
void getData(vector<vector<double>> &tab, vector<double> &b, int &taille);
void affiche(vector<vector<double>> tab);
void MethodeCholesky(vector<vector<double>>& A,vector<double>& B, int taille);
void verifierMatrice(vector<vector<double>>& A, int taille);
vector<vector<double>> transposeMatrice(vector<vector<double>>& A, int taille);
vector<double> ResolutionPremierEquation(vector<vector<double>>& A,vector<double>& b,int taille);
int main()
{
    vector<vector<double>> tab;
    vector<double> b;
    vector<double> x;
    int taille;
    getData(tab, b, taille); // on prend les données insérer dans le data.txt
    verifierMatrice(tab,taille);
    affichetab(tab, b);
    cout << endl;
    MethodeCholesky(tab,b,taille);
    x=ResolutionPremierEquation(tab,b,taille);
    solution(taille, tab, x);
    return 0;
}

void affichetab(vector<vector<double>> tab, vector<double> b)
{
    cout << "voici les équations à résoudre: \n" << endl;
    unsigned int i = 0, j = 0;
    for (i = 0; i < tab[0].size(); i++)
    {
        cout << "|\t";
        for (j = 0; j < tab[0].size(); j++)
        {
            cout << tab[i][j] << "\t";
        }
        cout << "| |\tX" << i + 1 << "\t|="
             << "|\t" << b[i] << "\t|";
        cout << endl;
    }
}

void echange(double *a, double *b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

void solution(int taille, vector<vector<double>> &tab, vector<double> &b)
{
    tab = transposeMatrice(tab,taille);
    // voici la fonction qui résout et affiche les résultats des équations
    double s(0);
    vector<double> x;
    for (int i = 0; i < taille; i++)
        x.push_back(0);
    for (int i = taille - 1; i >= 0; i--)
    {
        s = 0;
        for (int j = i + 1; j < taille; j++)
            s += (tab[i][j] * x[j]);
        x[i] = (b[i] - s) / tab[i][i];
    }
    cout << "voici les solutions aux équations: " << endl;
    for (int i = 0; i < taille; i++)
        cout << "X" << i + 1 << "=" << x[i] << endl;
}

void getData(vector<vector<double>> &tab, vector<double> &b, int &taille)
{
    // c'est une fonction qui lit les données entrer dans le fichiers.txt
    ifstream dossier("data.txt");
    if (dossier)
    {
        dossier >> taille; // on prend on premier lieu la taille de la matrice
        for (int j = 0; j < taille; j++)
        {
            vector<double> ligne;
            for (int i = 0; i < taille; i++)
            {
                double temp = 0;
                dossier >> temp;
                ligne.push_back(temp);
            }
            tab.push_back(ligne); // on insert les éléments de la ligne dans le tableau
        }
        for (int i = 0; i < taille; i++)
        {
            double temp = 0;
            dossier >> temp;
            b.push_back(temp); // on insert les éléments lignes
        }
    }
    else
    {
        cout << "les données à lire n'existe pas";
    }
}
void affiche(vector<vector<double>> tab)
{
    for (unsigned int i = 0; i < tab[0].size(); i++)
    {
        for (unsigned int j = 0; j < tab[0].size(); j++)
        {
            cout << "\t" << tab[i][j];
        }
        cout << endl;
    }
}

void MethodeCholesky(vector<vector<double>>& A,vector<double>& B, int taille){
	double sum = 0;
	for(int i(0); i<taille ; i++){
	    for(int j(0); j<taille ; j++){
				sum = 0;
				if(j>i){
					A[i][j]=0;
				}
				else if(i==j){
					for(int k(0); k<j ; k++){
						sum += A[i][k]*A[i][k];
					}
					A[i][i] = sqrt(A[i][i] - sum);
				
				}else{
					for(int k(0);k<j ; k++){
						sum += A[i][k]*A[j][k];
					}
					A[i][j] = 1.0 / A[j][j] * (A[i][j] - sum); 
		        }
			
		    }	
	    }
}


vector<double> ResolutionPremierEquation(vector<vector<double>>& A,vector<double>& b,int taille){
    double s(0);
    vector<double> x;
    for (int i = 0; i < taille; i++){
        x.push_back(0);
    }
        
    for (int i = 0; i < taille; i++)
    {
        s = 0;
        for (int j = 0; j < taille; j++)
            s += (A[i][j] * x[j]);
        x[i] = (b[i] - s) / A[i][i];
    }
    return x;
}

vector<vector<double>> transposeMatrice(vector<vector<double>>& A, int taille){
    vector<vector<double>> temp;
    vector<double>x;
    for(int i=0;i<taille;i++){
        x.clear();
        for(int j=0;j<taille;j++){
            x.push_back(A[j][i]);
        }
        temp.push_back(x);
    }
    return temp;
}

void verifierMatrice(vector<vector<double>>& A, int taille){
    vector<vector<double>> temp;
    temp = transposeMatrice(A,taille);
    for(int i= 0;i<taille;i++){
        for(int j=0;j<taille;j++){
            if(temp[i][j]-A[i][j]!=0){
                cout << "\n la matrice entrer n'est pas symétrique donc on ne peut pas utiliser la méthode de Cholesky"<< endl;
                cout << "\t essayer une autre matrice \n"<<endl;
                exit(0);
            }
        }
    }
}