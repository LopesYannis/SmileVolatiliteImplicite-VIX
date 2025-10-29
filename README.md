Avant d'exécuter,  remplacer en ligne 13 et 136 le lien vers la base de donnée "sp-index.txt" pour l'extraction.

La dernière fonction, GraphVolImpliciteSPX , permet de tracer les Smiles de volatilité implicite portant sur options sur indice SP500 en fonction de la maturité.
Cette fonction prend en compte un taux de dividende q=0.0217. 
Le prix de marché est déterminé comme étant la moyenne du Ask et Bid pour les calls et puts européens extraits de la base de données.
La maturité T, le strike K ainsi que le taux d'intérêt sans rique r sont égalements extraits de la base de donnée.
La formule de Black & Scholes portant sur les options européenne est utilisée pour trouver la volatilité implicite grâce à l'algorithme de Newton.

La fonction newtonTest trace un smile de volatilité avec des données de marché déclarées sous forme de liste dans la fonction.
