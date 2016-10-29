package algorithms;

import java.awt.Point;
import java.util.ArrayList;

public class Tree2D {
  private Point root;
  private ArrayList<Tree2D> subtrees;

  public Tree2D (Point p, ArrayList<Tree2D> trees){
    this.root=p;
    this.subtrees=trees;
  }
  public Point getRoot(){
    return this.root;
  }
  public ArrayList<Tree2D> getSubTrees(){
    return this.subtrees;
  }
  public double distanceRootToSubTrees(){
    double d=0;
    for (int i=0;i<this.subtrees.size();i++) d+=subtrees.get(i).getRoot().distance(root);
    return d;
  }
  
  
  public Tree2D Tree_clone() {
		if(this == null ) return null;
		ArrayList<Tree2D> liste = (ArrayList<Tree2D>) new ArrayList<Tree2D>(getSubTrees()).clone();

		return new Tree2D(getRoot(), liste);
	}
}
