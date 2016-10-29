package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;

public class Edge implements Comparable<Edge> {

	private Point p1;
	private Point p2;
	double distance;
	
	public Edge(Point p1, Point p2){
		this.p1=p1;
		this.p2=p2;
		distance=p1.distance(p2);
	}

	public Point getP1() {
		return p1;
	}


	public Point getP2() {
		return p2;
	}


	public double getDistance() {
		return distance;
	}



	@Override
	public int compareTo(Edge o) {
		// TODO Auto-generated method stub
		
		Double myDistance=distance;
		Double oDistance=o.getDistance();
		return myDistance.compareTo(oDistance) ;
	}
	
	
	
	/*int partition(ArrayList<Edge> edge,  int debut,  int fin)
	{
	   int compteur = debut;
	  double pivot =  edge.get(debut).getDistance();
	   int i;

	   for(i=debut+1; i<=fin; i++)
	   {
	      if(edge.get(i).getDistance() <pivot) // si élément inférieur au pivot
	      {
	         compteur++; // incrémente compteur cad la place finale du pivot
	         echanger(edge, compteur, i); // élément positionné
	      }
	   echanger(edge, compteur, debut); // le pivot est placé
	   return compteur; // et sa position est retournée
	}

	void triRapideAux(int tableau[], const int debut, const int fin)
	{
	   if(debut<fin) // cas d'arrêt pour la récursivité
	   {
	      int pivot = partition(tableau, debut, fin); // division du tableau
	      triRapideAux(tableau, debut, pivot-1); // trie partie1
	      triRapideAux(tableau, pivot+1, fin); // trie partie2
	   }
	}

	void triRapide(int tableau[], const int longueur)
	{
	   triRapideAux(tableau, 0, longueur-1);
	} 
	
}

	private void echanger(ArrayList<Edge> edge, int compteur, int i) {
		// TODO Auto-generated method stub
		
	}*/
}
