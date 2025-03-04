/*****************************************************************************
 * FILE NAME: monte_carlo_sphere.c                                           *
 *                                                                           *
 * COMPILATION:   gcc -Wall -o prog *.c -lm                                  *
 *                ./prog                                                     *
 *                                                                           *
 * PURPOSE:                                                                  *
 * Ce programme utilise la méthode de Monte Carlo pour estimer le            *
 * volume d'une sphère de rayon 1 et calcule un intervalle de                *
 * confiance basé sur plusieurs simulations.                                 *
 *                                                                           *
 * FILE REFERENCES:                                                          *
 * mt64.h - Contient la déclaration de la fonction genrand64_real1           *
 *          pour générer des nombres aléatoires.                             *
 *                                                                           *
 * EXTERNAL REFERENCES:                                                      *
 * genrand64_real1() - Fonction pour générer des nombres aléatoires.         *
 *                                                                           *
 * NOTES:                                                                    *
 * - La valeur critique t utilisée pour l'intervalle de confiance à          *
 *   95% est basée sur 19 degrés de liberté (pour 20 simulations).           *
 *                                                                           *
 * DEVELOPMENT HISTORY:                                                      *
 *    Date : 28/02/2025      Author: MOTA Titouan                            *
 *                                                                           *
 * ALGORITHM:                                                                *
 * 1. Estimer le volume de la sphère en utilisant la méthode de Monte Carlo. *
 * 2. Exécuter plusieurs simulations pour obtenir une moyenne des résultats. *
 * 3. Calculer l'intervalle de confiance basé sur la variance des résultats. *
 *                                                                           *
 *****************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mt64.h" 


/************************************************************
 *                                                          *
 * NOM DE LA FONCTION: estimateSphereVolumeMonteCarlo       *
 *                                                          *
 * ARGUMENTS:                                               *
 * ----------                                               *
 * numIterations : double                                   *
 *     Nombre d'itérations pour la simulation Monte Carlo.  *
 *                                                          *
 * RETOURNE:                                                *
 * ---------                                                *
 * double : Volume estimé de la sphère de rayon 1.          *
 *                                                          *
 ***********************************************************/

double estimateSphereVolumeMonteCarlo(double numIterations) {
    int         pointsInsideSphere = 0;           // Compteur pour les points à l'intérieur de la sphère
    double      randomX, randomY, randomZ;        // Coordonnées des points aléatoires
    double      distanceFromOrigin;               // Distance du point à l'origine

    // Boucle sur le nombre d'itérations
    for (int i = 0; i < numIterations; i++)
    {
        // Génération de points aléatoires
        randomX = genrand64_real1();
        randomY = genrand64_real1();
        randomZ = genrand64_real1();

        // Calcul de la distance à l'origine
        distanceFromOrigin = sqrt(pow(randomX, 2) + pow(randomY, 2) + pow(randomZ, 2));

        // Vérification si le point est à l'intérieur de la sphère de rayon 1
        if (distanceFromOrigin <= 1)
        {
            pointsInsideSphere += 1;
        }
    }

    // Retourne le volume estimé en ajustant par le facteur 8
    return ((double)pointsInsideSphere / numIterations) * 8;

    /*
    Pour avoir PI

    double volumeSphere = ((double)pointsInsideSphere / numIterations) * 8;
    return (volumeSphere * 3) / 4;
    */

}

/*********************************************************************
 *                                                                   *
 * NOM DE LA FONCTION: runMonteCarloSimulations                      *
 *                                                                   *
 * ARGUMENTS:                                                        *
 * ----------                                                        *
 * numPoints : double                                                *
 *     Nombre de points utilisés dans chaque simulation Monte Carlo. *
 *                                                                   *
 * numSimulations : double                                           *
 *     Nombre de simulations Monte Carlo à effectuer.                *
 *                                                                   *
 * RETOURNE:                                                         *
 * ---------                                                         *
 * double : Moyenne des résultats des simulations Monte Carlo.       *
 *                                                                   *
 ********************************************************************/

double runMonteCarloSimulations(double numPoints, double numSimulations){
    long double sumOfResults = 0;               // Somme des résultats des simulations

    // Boucle sur le nombre de simulations
    for (int i = 0; i < numSimulations; i++)
    {
        sumOfResults += estimateSphereVolumeMonteCarlo(numPoints);
    }

    // Retourne la moyenne des résultats des simulations
    return sumOfResults / numSimulations;
}

/**********************************************************************
 *                                                                    *
 * NOM DE LA FONCTION: calculateConfidenceInterval                    *
 *                                                                    *
 * ARGUMENTS:                                                         *
 * ----------                                                         *
 * numSimulations : double                                            *
 *     Nombre de simulations Monte Carlo à effectuer.                 *
 *                                                                    *
 * numPoints : double                                                 *
 *     Nombre de points utilisés dans chaque simulation Monte Carlo.  *
 *                                                                    *
 * RETOURNE:                                                          *
 * ---------                                                          *
 * Tableau contenant les bornes de l'intervalle de confiance          *
 *                                                                    *
 **********************************************************************/

double * calculateConfidenceInterval(double numSimulations, double numPoints){

    double * confidenceInterval = malloc(2 * sizeof(double));        // Allocation de mémoire pour stocker les bornes de l'intervalle de confiance
    double *results = malloc(numSimulations * sizeof(double));       // Allocation de mémoire pour stocker les résultats des simulations
    double   moyenne = 0;                                            // Variable pour stocker la moyenne des résultats des simulations
    double   varianceSum = 0;                                        // Variable pour accumuler les écarts par rapport à la moyenne

    //Verif allocation
    if(confidenceInterval == NULL || results == NULL)
    {
        fprintf(stderr, "Erreur d'allocation de mémoire\n");
        return NULL;
    }

    // Boucle pour accumuler les résultats de Monte Carlo
    for (int i = 0; i < numSimulations; i++)
    {   
        // Exécution de la simulation Monte Carlo
        results[i] = estimateSphereVolumeMonteCarlo(numPoints); 
        // Accumulation des résultats pour le calcul de la moyenne
        moyenne += results[i]; 
    }

    // Calcul de la moyenne
    moyenne = moyenne / numSimulations;

    // Calcul des écarts par rapport à la moyenne
    for (int i = 0; i < numSimulations; i++)
    {
        // Accumulation des écarts au carré
        varianceSum += (results[i] - moyenne) * (results[i] - moyenne); 
    }

    //Calcul de la variance
    double variance = varianceSum / (numSimulations - 1);

    //Calcul de l'intervalle de confiance avec 2.704 d'apres la table de la loi de student
    double confidenceRadius = 2.704 * sqrt(variance / numSimulations);

    confidenceInterval[0] = moyenne - confidenceRadius;
    confidenceInterval[1] = moyenne + confidenceRadius;

    printf("L'intervalle de confiance à 99%% est : [%.4f, %.4f]\n", confidenceInterval[0], confidenceInterval[1]);

    free(results);
    return confidenceInterval;
    
}


int main(){
   
    return 0;
}

    /* 
    
    Question  1
    printf("Volume estimé de la sphère (n itérations) : %f\n", estimateSphereVolumeMonteCarlo(1000000));
    
    Question 2 
    printf("Résultat moyen de n simulations (k points chacune) : %f\n", runMonteCarloSimulations(1000000, 40));

    Question 3
    double * confidenceInterval = calculateConfidenceInterval(40,1000000);
    free(confidenceInterval);
    
    */

