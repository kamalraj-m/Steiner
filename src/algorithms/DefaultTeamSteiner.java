package algorithms;

import java.awt.Point;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class DefaultTeamSteiner {

	
  public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
	  ArrayList<Point> hitPointFinal = new ArrayList<Point>();
	  ArrayList<Point> hitPointExtrem = new ArrayList<Point>();
	  ArrayList<Edge> edge = new ArrayList<Edge>();
	  ArrayList<Edge> kruskal =  kruskalEd(hitPoints, edge);
	  ArrayList<Point> hitClo = (ArrayList<Point>) hitPoints.clone();
	
	
		for(Point x : hitPoints){
		
			if(voisins(x, kruskal).size() == 1){
				hitPointExtrem.add(x);
			}
		}
	
		Point d = PMaxEx(hitPointExtrem, kruskal);
		for(Point p: hitPoints){
			if(!p.equals(d)) hitPointFinal.add(p);
		}
		
		hitPointExtrem.remove(d);
		
		Tree2D res =   calculSteiner(points, edgeThreshold, hitPointFinal);
	
		double scN = score(res);
		if(scN > 1664){
			res = calculSteinerBudget(points, edgeThreshold, hitPointFinal);
			scN = score(res);
		}
		
		return res;

  }
  
  
  public Point Point3OrdoneMin(ArrayList<Point> points){
		ArrayList<Integer> res = new ArrayList<Integer>();
		
	  for(int i = 0; i < points.size() ; i++){
		  res.add(points.get(i).y);
	  }
		  
	Collections.sort(res);
	Point pres = points.get(0);
	for(Point p: points){
		if(p.y == res.get(2)) pres = p;
	}
	
	
	return pres;
  }
  
  public Point DisMax(Point p, ArrayList<Point> points){
	  
	  double disMax = Double.MIN_VALUE;
	  Point res = points.get(0);
	  System.out.println("IZZIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
	  for(Point d: points){
		  if(d.distance(p) > disMax ){
			  disMax = d.distance(p);
			  res = d;
					  
		  }
	  }
	  
	  return res;
  }
  
  public Point PMaxEx(ArrayList<Point> p, ArrayList<Edge> edge){
	 double disMax = Double.MIN_VALUE;
	 Point res = p.get(0);
	for(Point d: p){
		if((voisins(d, edge).get(0).distance(d) >disMax)){
			disMax = voisins(d, edge).get(0).distance(d);
			res = d;
		}
	}
	
	return res;
  }
  

	public int[][] floydWarshall(ArrayList<Point> points, int threshhold) {
		Double[][] D = new Double[points.size()][points.size()];
		int[][] P = new int[points.size()][points.size()];

		for (int i = 0; i < points.size(); i++)
			for (int j = 0; j < points.size(); j++) {
				Point p = points.get(i);
				Point q = points.get(j);
				double distance = p.distance(q);
				if (distance > threshhold) {
					D[i][j] = Double.POSITIVE_INFINITY;
				} else {
					D[i][j] = distance;
				}
				P[i][j] = j;
			}

		for (int k = 0; k < D.length; k++) {
			for (int i = 0; i < D.length; i++) {
				for (int j = 0; j < D.length; j++) {
					if (D[i][k] + D[k][j] < D[i][j]) {
						D[i][j] = D[i][k] + D[k][j];
						P[i][j] = P[i][k];
					}
				}
			}
		}

		return P;
	}

	double minDist(double d1, double d2) {

		return Math.min(d1, d2);
	}

	public int[][] calculShortestPaths(ArrayList<Point> points,
			int edgeThreshold) {
		double[][] D = new double[points.size()][points.size()];
		int[][] P = new int[points.size()][points.size()];

		for (int i = 0; i < points.size(); i++)
			for (int j = 0; j < points.size(); j++) {
				Point p = points.get(i);
				Point q = points.get(j);
				double distance = p.distance(q);
				if (distance > edgeThreshold) {
					D[i][j] = Double.POSITIVE_INFINITY;
				} else {
					D[i][j] = distance;
				}
				P[i][j] = j;
			}

		for (int k = 0; k < D.length; k++) {
			for (int i = 0; i < D.length; i++) {
				for (int j = 0; j < D.length; j++) {
					if (D[i][k] + D[k][j] < D[i][j]) {
						D[i][j] = D[i][k] + D[k][j];
						P[i][j] = P[i][k];
					}
				}
			}
		}
		return P;
	}
	
	public Point distMin(ArrayList<Point> points,Point z){
		Point res = points.get(0);
		double disMin = Double.MAX_VALUE;
		for(Point p: points){
			if(p.equals(z)) {return p;}
			else if(p.distance(z) < disMin){
				disMin = p.distance(z);
				res = p;
			}
		}
		
		return res;
	}

	public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold,
			ArrayList<Point> hitPoints) {
		
	
	return calculSteinerTempOp(points, edgeThreshold, hitPoints);
		
	}
	
	public Tree2D calculSteinerNew(ArrayList<Point> points, int edgeThreshold,
			ArrayList<Point> hitPoints, Tree2D ancien){
		
		ArrayList<Edge> ledge1 = new ArrayList<Edge>();
		ArrayList<Edge> kruskal =  kruskalEd(hitPoints, ledge1);
		
		for(Point p: hitPoints){
			ArrayList<Point> voisins = voisins(p, kruskal);

			if(voisins.size() < 3){
				continue;
			}
			else{
				for(int i= 0; i< voisins.size()-1; i++){
					for(int j= i; j< voisins.size(); j++){
						if(i == j) continue;
						for(int k= i; k< voisins.size(); k++){
						if(k == j) continue;
						System.out.println("j===>"+j);
						Point QPoint = Avec4Point(p, voisins.get(i), voisins.get(j) , voisins.get(k));
						for(int l = 0 ; l<15; l++){
							Point coté = distMin(points, QPoint);

						hitPoints.add(coté);
						Tree2D tmp =  calculSteinerIni(points, edgeThreshold, hitPoints); 
						double newscore = score(tmp);
						if(newscore > score(ancien)){
							System.out.println("nouveau score " +newscore);
							QPoint = coté;
							hitPoints.remove(coté);
						}
						else{
							System.out.println("HitPointsTemp aprés l'ajout == "+hitPoints.size());
							System.out.println("new Score == "+ score(tmp));
					
							return tmp;
						}

						}
					}
					}
				}
			}
		}
		return ancien;
	}
	
	
	
	public Tree2D calculSteinerIni(ArrayList<Point> points, int edgeThreshold,
			ArrayList<Point> hitPoints){
	ArrayList<Edge> ledge = new ArrayList<Edge>();
	
	ArrayList<Edge> kruskal =  kruskalEd(hitPoints, ledge); // la liste des arretes aprés kruskal
	double[][] D = new double[points.size()][points.size()];
	int[][] P = new int[points.size()][points.size()];

	for (int i = 0; i < points.size(); i++)
		for (int j = 0; j < points.size(); j++) {
			Point p = points.get(i);
			Point q = points.get(j);
			double distance = p.distance(q);
			if (distance > edgeThreshold) {
				D[i][j] = Double.POSITIVE_INFINITY;
			} else {
				D[i][j] = distance;
			}
			P[i][j] = j;
		}

	for (int k = 0; k < D.length; k++) {
		for (int i = 0; i < D.length; i++) {
			for (int j = 0; j < D.length; j++) {
				if (D[i][k] + D[k][j] < D[i][j]) {
					D[i][j] = D[i][k] + D[k][j];
					P[i][j] = P[i][k];
				}
			}
		}
	}
	
	ArrayList<Edge> aretes = new ArrayList<Edge>();
	ArrayList<Point> puntos = new ArrayList<Point>();
	for (int i = 0; i < kruskal.size(); i++) {
		int p1 = points.indexOf(kruskal.get(i).getP1());
		int p2 = points.indexOf(kruskal.get(i).getP2());
		ArrayList<Edge> tmpEdge = new ArrayList<Edge>();
		ArrayList<Point> tmpPoints = new ArrayList<Point>();
		while(p1 != p2 ) {
			int k = P[p1][p2];
			
			
			tmpEdge.add(new Edge(points.get(p1), points.get(k)));
			tmpPoints.add(points.get(k));
			
			p1 = k;
		}
		aretes.addAll(tmpEdge);
		puntos.addAll(tmpPoints);
		
	
	}
	hitPoints.addAll(puntos);
	Tree2D res =  kruskal(hitPoints , aretes);
	return res;
}
	
	
	
	
	public Tree2D calculSteinerTempOp(ArrayList<Point> points, int edgeThreshold,
			ArrayList<Point> hitPoints) {
		
		 
		 Tree2D  temp = recher_local(kruskalNew(hitPoints), new ArrayList<Point>(),points);		 double[][] D = new double[points.size()][points.size()];
			int[][] P = new int[points.size()][points.size()];

			for (int i = 0; i < points.size(); i++)
				for (int j = 0; j < points.size(); j++) {
					Point p = points.get(i);
					Point q = points.get(j);
					double distance = p.distance(q);
					if (distance > edgeThreshold) {
						D[i][j] = Double.POSITIVE_INFINITY;
					} else {
						D[i][j] = distance;
					}
					P[i][j] = j;
				}

			for (int k = 0; k < D.length; k++) {
				for (int i = 0; i < D.length; i++) {
					for (int j = 0; j < D.length; j++) {
						if (D[i][k] + D[k][j] < D[i][j]) {
							D[i][j] = D[i][k] + D[k][j];
							P[i][j] = P[i][k];
						}
					}
				}
			}
			 double newScore = 0;
			 Tree2D min = null;	
			/*rajouter les points bleu*/
			ArrayList<Point> l = new ArrayList<Point>(); 
			ArrayList<Point> res = new ArrayList<Point>();
			res.add(temp.getRoot());

			if(temp.getSubTrees().size()==0){
				l = res;
			}

			for (Tree2D  tree : temp.getSubTrees()){
				Point p = temp.getRoot();
				Point q = tree.getRoot();

				Point K = p;
				int pI = points.indexOf(p);
				int qI = points.indexOf(q);

				while(! K.equals(q)){

					K = points.get(P[pI][qI]);
					if(! K.equals(q))
						res.add(K);
					pI = points.indexOf(K);

				}

				res.addAll(AjoutBluePoint(points, P, tree));
			}

		 l = res;
		 /*fin rajoue de point bleu*/
		 
		newScore = score(recher_local(kruskalNew(l), new ArrayList<Point>(), points));
			min = recher_local(kruskalNew(l), new ArrayList<Point>(), points);
		


		 
		 return  min;

		
	}
	
	
	
	public ArrayList<Point> BluePointsToAdd(ArrayList<Point> points, int [][] path, Tree2D redTree) {

		ArrayList<Point> res = new ArrayList<Point>();
		res.add(redTree.getRoot());

		if(redTree.getSubTrees().size()==0){
			return res;
		}

		for (Tree2D  tree : redTree.getSubTrees()){
			Point p = redTree.getRoot();
			Point q = tree.getRoot();

			Point TMP = p;
			int pI = points.indexOf(p);
			int qI = points.indexOf(q);

			while(! TMP.equals(q)){

				TMP = points.get(path[pI][qI]);
				if(! TMP.equals(q))
					res.add(TMP);
				pI = points.indexOf(TMP);

			}

			res.addAll(BluePointsToAdd(points, path, tree));
		}


		return res;
	}
	
	//Compte score 
	
	
	
	public double score(Tree2D t){
		double resultat = 0;
		if(t.getSubTrees() != null){
			for (int i = 0; i< t.getSubTrees().size(); i++){
				resultat +=score(t.getSubTrees().get(i));
		}
			return t.distanceRootToSubTrees() + resultat;
		}
		else{
			return 0;
		}
	}
	
	
	
	//Baricentre 
		public Point bary(Point a, Point b, Point c){
			int xG = (a.x + b.x + c.x)/3;
			int yG = (a.y + b.y + c.y)/3;
			return new Point (xG, yG);
			
		}
		
		//Point cà coté	
	public ArrayList<Point> voisins(Point p, ArrayList<Edge> edge){
		ArrayList<Point> resultat = new ArrayList<Point>();
		
		for(Edge e: edge){
			if( e.getP1() == p){
				resultat.add(e.getP2());
			}
			else if(e.getP2() == p){
				resultat.add(e.getP1());
			}
			else continue;
		}
		
		return resultat;
	}

	
/**
* Kruskal avec un retour de Tree2D
* @param pts
* @return
*/
		public Tree2D kruskal(ArrayList<Point> pts ,ArrayList<Edge> ledge){
			// ArrayList<Point> points = new ArrayList<Point>(pts);
		//	ArrayList<Edge> ledge = new ArrayList<Edge>();
			ArrayList<Point> res = new ArrayList<Point>(pts);
			HashMap<Point, Integer> pet = new HashMap<Point, Integer>();
			
			int compteur = 0;
			for (Point p : pts)
				pet.put(p, -1);

			ArrayList<Edge> tmp = new ArrayList<Edge>();
			for (Point p : pts)
				for (Point q : pts) {
					if (p.equals(q))
						continue;
					tmp.add(new Edge(p, q));
				}

			Collections.sort((List<Edge>) tmp);

	ledge = (ArrayList<Edge>) tmp.clone();
			ArrayList<Edge> solutions = new ArrayList<Edge>();
			// solutions.add(ledge.get(0));
			int j = 0;
			for (Edge e : ledge) {

				// cycleExists(solutions, e, pet,pts);
				if (pet.get(e.getP1()) == -1 && pet.get(e.getP2()) == -1) {
					
					pet.put(e.getP1(), compteur);
					pet.put(e.getP2(), compteur);
					solutions.add(e);
					compteur++;

				} else if (pet.get(e.getP1()) == -1 && pet.get(e.getP2()) != -1) {

					pet.put(e.getP1(), pet.get(e.getP2()));
					solutions.add(e);

				} else if (pet.get(e.getP1()) != -1 && pet.get(e.getP2()) == -1) {

					pet.put(e.getP2(), pet.get(e.getP1()));
					solutions.add(e);

				} else if (pet.get(e.getP1()) == pet.get(e.getP2())) {

				} else {
					if (pet.get(e.getP1()) != pet.get(e.getP2())) {

						for (Point p : pts) {
							if (p.equals(e.getP1()))
								continue;
							if (pet.get(e.getP1()) == pet.get(p))
								pet.put(p, pet.get(e.getP2()));
						}
						pet.put(e.getP1(), pet.get(e.getP2()));

						solutions.add(e);
					}
				}

			}
			Point first = solutions.get(0).getP1();
			System.out.println(solutions.size());
			return convert(first, solutions);

		}
	

		
		
		
		
/**
* Krustal retour liste de edge		
* @param points
* @return
*/
		
		public ArrayList<Edge> kruskalEd(ArrayList<Point> points , ArrayList<Edge> ledge){
			// ArrayList<Point> points = new ArrayList<Point>(pts);
		//	ArrayList<Edge> ledge = new ArrayList<Edge>();
			ArrayList<Point> res = new ArrayList<Point>(points);
			HashMap<Point, Integer> pet = new HashMap<Point, Integer>();
			
			int compteur = 0;
			for (Point p : points)
				pet.put(p, -1);

			ArrayList<Edge> tmp = new ArrayList<Edge>();
			for (Point p : points)
				for (Point q : points) {
					if (p.equals(q))
						continue;
					tmp.add(new Edge(p, q));
				}

			Collections.sort((List<Edge>) tmp);

	ledge = (ArrayList<Edge>) tmp.clone();

			ArrayList<Edge> solutions = new ArrayList<Edge>();
			int j = 0;
			for (Edge e : ledge) {
				if (pet.get(e.getP1()) == -1 && pet.get(e.getP2()) == -1) {
					
					pet.put(e.getP1(), compteur);
					pet.put(e.getP2(), compteur);
					solutions.add(e);
					compteur++;

				} else if (pet.get(e.getP1()) == -1 && pet.get(e.getP2()) != -1) {

					pet.put(e.getP1(), pet.get(e.getP2()));
					solutions.add(e);

				} else if (pet.get(e.getP1()) != -1 && pet.get(e.getP2()) == -1) {

					pet.put(e.getP2(), pet.get(e.getP1()));
					solutions.add(e);

				} else if (pet.get(e.getP1()) == pet.get(e.getP2())) {

				} else {
					if (pet.get(e.getP1()) != pet.get(e.getP2())) {

						for (Point p : points) {
							if (p.equals(e.getP1()))
								continue;
							if (pet.get(e.getP1()) == pet.get(p))
								pet.put(p, pet.get(e.getP2()));
						}
						pet.put(e.getP1(), pet.get(e.getP2()));

						solutions.add(e);
					}
				}

			}
			Point first = solutions.get(0).getP1();
			return solutions;

		}
		
		
	
		
		public Tree2D kruskalNew(ArrayList<Point> points){
			// ArrayList<Point> points = new ArrayList<Point>(pts);
			ArrayList<Edge> ledge = new ArrayList<Edge>();
			ArrayList<Point> res = new ArrayList<Point>(points);
			HashMap<Point, Integer> pet = new HashMap<Point, Integer>();
			
			int compteur = 0;
			for (Point p : points)
				pet.put(p, -1);

			ArrayList<Edge> tmp = new ArrayList<Edge>();
			for (Point p : points)
				for (Point q : points) {
					if (p.equals(q))
						continue;
					ledge.add(new Edge(p, q));
				}

			Collections.sort((List<Edge>) ledge);

			ArrayList<Edge> solutions = new ArrayList<Edge>();
			int j = 0;
			for (Edge e : ledge) {
				if (pet.get(e.getP1()) == -1 && pet.get(e.getP2()) == -1) {
					
					pet.put(e.getP1(), compteur);
					pet.put(e.getP2(), compteur);
					solutions.add(e);
					compteur++;

				} else if (pet.get(e.getP1()) == -1 && pet.get(e.getP2()) != -1) {

					pet.put(e.getP1(), pet.get(e.getP2()));
					solutions.add(e);

				} else if (pet.get(e.getP1()) != -1 && pet.get(e.getP2()) == -1) {

					pet.put(e.getP2(), pet.get(e.getP1()));
					solutions.add(e);

				} else if (pet.get(e.getP1()) == pet.get(e.getP2())) {

				} else {
					if (pet.get(e.getP1()) != pet.get(e.getP2())) {

						for (Point p : points) {
							if (p.equals(e.getP1()))
								continue;
							if (pet.get(e.getP1()) == pet.get(p))
								pet.put(p, pet.get(e.getP2()));
						}
						pet.put(e.getP1(), pet.get(e.getP2()));

						solutions.add(e);
					}
				}

			}
			Point first = solutions.get(0).getP1();
			System.out.println(solutions.size());
			return convert(first, solutions);

		}
		
		
		public ArrayList<Edge> split(ArrayList<Edge> solution, ArrayList<Edge> l) {
			ArrayList<Edge> tmp = (ArrayList<Edge>) solution.clone();
			for (Edge e : l) {
				if (tmp.contains(e))
					tmp.remove(e);
			}
			return tmp;
		}

		public Tree2D convert(Point p, ArrayList<Edge> solution) {
			if (solution.size() == 0) {
				return new Tree2D(p, new ArrayList<Tree2D>());
			} else {
				ArrayList<Tree2D> fils = new ArrayList<Tree2D>();

				ArrayList<Edge> l = new ArrayList<Edge>();
				for (Edge e : solution) {
					if (e.getP1().equals(p) || e.getP2().equals(p))
						l.add(e);
				}

				for (Edge e : l) {
					if (e.getP1() == p)
						fils.add(convert(e.getP2(), split(solution, l)));
					if (e.getP2() == p) {
						fils.add(convert(e.getP1(), split(solution, l)));
					}
				}
				return new Tree2D(p, fils);
			}

		}
	 
		
		
/*
* 
* 
* NOUVEAU
* 		
*/

public  Tree2D recher_local(Tree2D arbre, ArrayList<Point> cigle, ArrayList<Point> points){

	ArrayList<Point> tmp = new ArrayList<Point>();


	if (!cigle.contains(arbre.getRoot())) {

		cigle.add(arbre.getRoot());
		tmp.add(arbre.getRoot());

		for (Tree2D x : arbre.getSubTrees()) {
			recher_local(x, cigle, points);
			tmp.add(x.getRoot());


		}

		if (arbre.getSubTrees().size() >1 ) {

			Point fer = PointFermat(tmp.get(0), tmp.get(1), tmp.get(2));
			Point p = P_distMin(fer, points); //Retourne le point de distance min avec le oint fermat

			if (D_LIST(p, tmp) < D_TREE(arbre)) { 
				Tree2D TreeFermat = new Tree2D(p, (ArrayList<Tree2D>) arbre.getSubTrees().clone());
				arbre.getSubTrees().clear();
				arbre.getSubTrees().add(TreeFermat);
			}
		}

		if (arbre.getSubTrees().size() >= 1) {

			for (int i=0; i < arbre.getSubTrees().size(); i++) {
				Tree2D t = arbre.getSubTrees().get(i);

				for (int j = 0; j < t.getSubTrees().size(); j++) {
					Tree2D tt = t.getSubTrees().get(j);

					tmp = new ArrayList<Point>();
					tmp.add(arbre.getRoot());
					tmp.add(t.getRoot());
					tmp.add(tt.getRoot());


					Point fer = PointFermat(tmp.get(0), tmp.get(1), tmp.get(2));
					Point p = P_distMin(fer, points);

					if (D_LIST(p, tmp) < D_LIST(tmp.get(1), tmp)) {
						ArrayList<Tree2D> listeTMP = new ArrayList<Tree2D>();
						listeTMP.add(t.Tree_clone());
						listeTMP.add(tt.Tree_clone());
						listeTMP.get(0).getSubTrees().remove(j);
						Tree2D TreeFermat = new Tree2D(p,(ArrayList<Tree2D>) listeTMP.clone());
						arbre.getSubTrees().remove(i);
						arbre.getSubTrees().add(i,TreeFermat);
					}
				}
			}
		}

	}
	return arbre;
}




public double D_TREE(Tree2D arbre){

	if(arbre == null || arbre.getSubTrees().size()==0) return 0;
	double res = 0;
	for (Tree2D x : arbre.getSubTrees()) {
		res += arbre.getRoot().distance(x.getRoot());
	}
	return res;
}

public double D_LIST(Point F, ArrayList<Point> liste){

	if(liste == null || liste.size()==0) return 0;
	double res = 0;
	for (Point p : liste) {
		res += F.distance(p);
	}
	return res;
}



public Point PointFermat(Point A, Point B, Point C) {


	double AB = A.distance(B);
	double BC = C.distance(B);
	double CA = A.distance(C);

	double angleA = angle(A,  B, C); 
	double angleB = angle(B,  A, C); 
	double angleC = angle(C,  B, A); 


	if( ((angleA * 180 / Math.PI) < 120) && ((angleB * 180 / Math.PI) < 120) && ((angleC * 180 / Math.PI) < 120)){

		Point F = new Point (     

				(int) ((BC * A.x / Math.sin(angleA+Math.PI/3) + 
						CA * B.x / Math.sin(angleB+Math.PI/3) +
						AB * C.x / Math.sin(angleC+Math.PI/3) 
						)/ (
								BC / Math.sin(angleA+Math.PI/3) +
								CA / Math.sin(angleB+Math.PI/3) +
								AB / Math.sin(angleC+Math.PI/3) 
								))

				,

				(int) ((BC * A.y / Math.sin(angleA+Math.PI/3) + 
						CA * B.y / Math.sin(angleB+Math.PI/3) +
						AB * C.y / Math.sin(angleC+Math.PI/3) 
						)/ (
								BC / Math.sin(angleA+Math.PI/3) +
								CA / Math.sin(angleB+Math.PI/3) +
								AB / Math.sin(angleC+Math.PI/3) 
								))





				);

		return F;
	}
	else if((angleA * 180 / Math.PI) >= 120) return A;
	else if((angleB * 180 / Math.PI) >= 120) return B;
	else return C;

}



public Point P_distMin(Point point, ArrayList<Point> points) {

	if(points.isEmpty())return null;
	Point result = points.get(0);
	for(Point p : points) if(p.distance(point)<=result.distance(point)) result = p;
	return result;

}


public double angle(Point C, Point B, Point A) {
	double a = C.distance(B);
	double b = C.distance(A);
	double c = B.distance(A);

	return Math.acos((Math.pow(a, 2) + Math.pow(b, 2) - Math.pow(c, 2))
			/ (2 * a * b));
}


public ArrayList<Point> AjoutBluePoint(ArrayList<Point> points, int [][] path, Tree2D redTree) {

	ArrayList<Point> res = new ArrayList<Point>();
	res.add(redTree.getRoot());

	if(redTree.getSubTrees().size()==0){
		return res;
	}

	for (Tree2D  tree : redTree.getSubTrees()){
		Point p = redTree.getRoot();
		Point q = tree.getRoot();

		Point TMP = p;
		int pI = points.indexOf(p);
		int qI = points.indexOf(q);

		while(! TMP.equals(q)){

			TMP = points.get(path[pI][qI]);
			if(! TMP.equals(q))
				res.add(TMP);
			pI = points.indexOf(TMP);

		}

		res.addAll(AjoutBluePoint(points, path, tree));
	}


	return res;
}


public Point Avec4Point(Point a, Point b, Point c, Point d ){
	int xG = (a.x + b.x + c.x + d.x)/4;
	int yG = (a.y + b.y + c.y + d.y)/4      ;
	return new Point (xG, yG);
	
}



public  Tree2D recher_local2(Tree2D arbre, ArrayList<Point> cigle, ArrayList<Point> points){

	ArrayList<Point> tmp = new ArrayList<Point>();


	if (!cigle.contains(arbre.getRoot())) {

		cigle.add(arbre.getRoot());
		tmp.add(arbre.getRoot());

		for (Tree2D x : arbre.getSubTrees()) {
			recher_local(x, cigle, points);
			tmp.add(x.getRoot());


		}

		if (arbre.getSubTrees().size() >1 ) {

			Point fer = PointFermat(tmp.get(0), tmp.get(1), tmp.get(2));
			Point p = P_distMin(fer, points); //Retourne le point de distance min avec le oint fermat

			if (D_LIST(p, tmp) < D_TREE(arbre)) { 
				Tree2D TreeFermat = new Tree2D(p, (ArrayList<Tree2D>) arbre.getSubTrees().clone());
				arbre.getSubTrees().clear();
				arbre.getSubTrees().add(TreeFermat);
			}
		}

		if (arbre.getSubTrees().size() >= 1) {

			for (int i=0; i < arbre.getSubTrees().size(); i++) {
				Tree2D t = arbre.getSubTrees().get(i);

				for (int j = 0; j < t.getSubTrees().size(); j++) {
					Tree2D tt = t.getSubTrees().get(j);

					tmp = new ArrayList<Point>();
					tmp.add(arbre.getRoot());
					tmp.add(t.getRoot());
					tmp.add(tt.getRoot());


					Point fer = Avec4Point(tmp.get(0), tmp.get(1), tmp.get(2),tmp.get(3));
					Point p = P_distMin(fer, points);

					if (D_LIST(p, tmp) < D_LIST(tmp.get(1), tmp)) {
						ArrayList<Tree2D> listeTMP = new ArrayList<Tree2D>();
						listeTMP.add(t.Tree_clone());
						listeTMP.add(tt.Tree_clone());
						listeTMP.get(0).getSubTrees().remove(j);
						Tree2D TreeFermat = new Tree2D(p,(ArrayList<Tree2D>) listeTMP.clone());
						arbre.getSubTrees().remove(i);
						arbre.getSubTrees().add(i,TreeFermat);
					}
				}
			}
		}

	}
	return arbre;
}
}
