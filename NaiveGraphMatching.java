import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeSet;

import org.neo4j.graphdb.Direction;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;

class NodeMatch {
	/**
	 * In the process of guessing potential nodes in solutions
	 */
	Node nodeMatched;
	int index;

	NodeMatch(Node n, int i) {
		nodeMatched = n;
		index = i;
	}

}


/***
 * finds the naive graph matched for the given target and query graph 
 * @author Lakshmi Ravi
 *
 */
public class NaiveGraphMatching {
	// private static String DB_PATH =
	// "C:\\Users\\laksh\\Documents\\Graph_Databases\\Neo4j\\default.graph.db\\";
	private static String DB_PATH;
	private GraphDatabaseService graphDb;

	// the query is denoted by the following data structures
	private List<String> queryNodes = new LinkedList<String>();
	// maps the name given to the node to its' label
	private HashMap<String, String> nodeLabelMap = new HashMap<String, String>();
	// maps the name given to a node to its' index - access
	private HashMap<String, Integer> nodeBucketIdMap = new HashMap<String, Integer>();

	// maps the name given to a node to it's profile Array
	private HashMap<Integer, List<String>> nodeProfileMap = new HashMap<Integer, List<String>>();

	private List<Node>[] nodeBuckets = null;
	private Integer[][] incomingEdges = null;
	private List<String[]> edgeLines = new LinkedList<String[]>();

	public ResourceIterator<Node>[] getNodes() {
		/**
		 * This is the method which defines the search space of the problem
		 * getNodes method filters the nodes of our interest
		 */
		ResourceIterator<Node>[] allNodes = new ResourceIterator[queryNodes.size()];
		int nodeIndex = 0;
		// get all the nodes in each of the bucket specified
		for (String node_name : queryNodes) {
			String labelForNode = nodeLabelMap.get(node_name);

			nodeBucketIdMap.put(node_name, nodeIndex);
			//create an empty profile when node-name is being Mapped to an index
			nodeProfileMap.put(nodeIndex, new LinkedList<String>());
			//append the current label of the node in it's profile 
			nodeProfileMap.get(nodeIndex).add(labelForNode);
			ResourceIterator<Node> nodesMatchingLabel = graphDb.findNodes(Label.label(labelForNode));
			allNodes[nodeIndex] = nodesMatchingLabel;
			nodeIndex++;
		}
		return allNodes;
		
	}
		
	
	
	/**
	 * From the bucket of all nodes that has matched the label - filter the nodes that has 
	 * @param allNodesMatchingLabel
	 * @return
	 */
		
		public List<Node>[] filterProfileBasedNodes(ResourceIterator<Node>[] allNodesMatchingLabel) {
			System.out.println("Filtering based on profile");
			List<Node>[] nodeBuckets = null;
			nodeBuckets = new LinkedList[queryNodes.size()];
			for(int i=0;i<queryNodes.size();i++){
				nodeBuckets[i] = new LinkedList<Node>();
			}
			//System.out.println("Filtering nodes based on profile for ..."+queryNodes.size()+" nodes");

			//iterate through all the nodes with same label and filter the one with profile matched
			int nodeIndex=0;
		for(ResourceIterator<Node> nodesMatchingLabel : allNodesMatchingLabel){
			
			//for all the nodes in that bucket
			while (nodesMatchingLabel.hasNext()) {
				Node n = nodesMatchingLabel.next();
				//get the profile of the current Node
				List<String> profileOfNode = nodeProfileMap.get(nodeIndex);
				//System.out.println(profileOfNode+" is the profile of"+nodeIndex);
				String[] profileArr = (String[]) n.getProperty("PROFILE");

				// if the profile of the current Node-'n' contains all the values as
				// in the profile of query node
				if (Arrays.asList(profileArr).containsAll(profileOfNode)) {
					// profile of the node matched- hence add to the selected bucket
					nodeBuckets[nodeIndex].add(n);
				}

			}
			nodeIndex++;
		}
		
		return nodeBuckets;
	}

	/***
	 * Utility : from the raw console input - create an adjacency list to
	 * represent edge relations between the nodes. to create a profile value for
	 * a node.
	 */
	private void setIncomingEdges() {
		Iterator<String[]> edge = edgeLines.iterator();
		int numberNodes = queryNodes.size();
		incomingEdges = new Integer[numberNodes][numberNodes];
		for (int i = 0; i < numberNodes; i++) {
			incomingEdges[i] = new Integer[numberNodes];
			for(int j=0;j<numberNodes;j++){
				incomingEdges[i][j] = new Integer(0);
			}
		}

		while (edge.hasNext()) {
			String[] oneEdge = edge.next();
			String start = oneEdge[0];
			String end = oneEdge[1];
			// to every outgoing edge - neighbor of a node, add the profile
			String neighborLabel = nodeLabelMap.get(end);

			int startIndex = this.nodeBucketIdMap.get(start);
			nodeProfileMap.get(startIndex).add(neighborLabel);

			int endIndex = this.nodeBucketIdMap.get(end);
			incomingEdges[startIndex][endIndex]=1;

		}

	}

	/**
	 * gets the nodes and edge list from the raw input lines of igraph file =>
	 * Query representation
	 * 
	 * @param inputLines
	 * @throws IOException
	 */
	public void getNodesEdgesIgraphQuery(File igraph) throws IOException {
		System.out.println("Getting Nodes and Edges from the IGraph file..."+igraph.getName());
		BufferedReader fileContent = new BufferedReader(new FileReader(igraph));
		String currentLine;

		while ((currentLine = fileContent.readLine()) != null) {
			String[] lineElements = currentLine.split(" ");
			if (lineElements[0].equals("e")) {
				// add a relation from the mentioned nodes
				String[] edgeLine = new String[2];
				edgeLine[0] = lineElements[1];
				edgeLine[1] = lineElements[2];
				edgeLines.add(edgeLine);

			} else if (lineElements[0].equals("v")) {
				queryNodes.add(lineElements[1]);
				nodeLabelMap.put(lineElements[1], lineElements[2]);
			}

		}
		fileContent.close();

	}

	/**
	 * gets the nodes and edge list from the raw input lines of igraph file =>
	 * Query representation
	 * 
	 * @param inputLines
	 * @throws IOException
	 */
	public void getNodesEdgesProteinQuery(File protein) throws IOException {
		System.out.println("Get nodes and edges from Protein file"+protein.getName());
		BufferedReader fileContent = new BufferedReader(new FileReader(protein));
		String currentLine;
		int nodeCount = Integer.parseInt(fileContent.readLine());

		// read the number of lines equal to nodeCount and add nodes in query
		while ((currentLine = fileContent.readLine()) != null
				&& nodeCount-- > 0) {
			String[] lineElements = currentLine.split(" ");
			queryNodes.add(lineElements[0]);
			nodeLabelMap.put(lineElements[0], lineElements[1]);

		}

		// for the rest of the lines add edges
		while ((currentLine = fileContent.readLine()) != null) {
			String[] lineElements = currentLine.split(" ");
			if(lineElements.length ==2){
			String[] edgeLine = new String[2];
			edgeLine[0] = lineElements[0];
			edgeLine[1] = lineElements[1];
			edgeLines.add(edgeLine);
			}

		}

		fileContent.close();
	}

	public void printSolution(HashMap<Integer, Node> solution) {
		System.out.println("\n***************************************");
		Iterator<Integer> solutionIndices = solution.keySet().iterator();
		while (solutionIndices.hasNext()) {
			int nodeIndex=  solutionIndices.next();
			String Nodename="";
			
			//iterate over the nodeBucketMap to map the keys 
			for (String name : nodeBucketIdMap.keySet()) {
			      if (nodeBucketIdMap.get(name).equals(nodeIndex)) {
			    	  	Nodename = name;
			      }
			    }
			Node sol_node = solution.get(nodeIndex);
			System.out.println(Nodename+":" + sol_node.getId());
		}
		System.out.println("***************************************");
	}

	public void startGraphMatching(int index, HashMap<Integer, Node> solution, int[] processingOrder) {
		// start the matching beginning from the node_bucket with the given
		// nodex
		//System.out.println("Processing Node index ==>"+index);
		int nodeIndexToProcess=-1;
		
		//if the processing order is not specified then 
		if (processingOrder !=null){
				nodeIndexToProcess = processingOrder[index];
		}
		else{
			nodeIndexToProcess = index;
		}
		List<Node> potential_i = this.nodeBuckets[nodeIndexToProcess];
		Iterator<Node> iNodeIterator = potential_i.iterator();

		// iterate through the edges of the query node
		Integer[]edges_i = this.incomingEdges[nodeIndexToProcess];
		List edges = new LinkedList<Integer>();
		for(int i=0;i<edges_i.length;i++){
			if(edges_i[i] ==1){
				edges.add(i);
			}
		}
		Iterator<Integer> edgeIterator = edges.iterator();
		int total_nodes_to_check = this.queryNodes.size();

		if (index >= total_nodes_to_check) {
			printSolution(solution);
			return;
		}

		while (iNodeIterator.hasNext()) {
			// for every edge listed in that find the sequence of nodes.
			Node currentNode = (Node) iNodeIterator.next();

			// check if the current-node potential_i matches with the solution
			// computed so-far

			// check if the node can be said as the 'index' and it doesn't
			// contradict with the solution computed so far
			if ( (!solution.containsValue(currentNode)) && (canMap(currentNode, nodeIndexToProcess, Arrays.asList(edges_i), solution))) {
				solution.put(nodeIndexToProcess, currentNode);
				startGraphMatching(index + 1, solution, processingOrder);
				solution.remove(currentNode);
			}
		}

	}

	/***
	 * currentNode can Map as the "index"
	 * 
	 * @param currentNode
	 * @param index
	 * @param edges_i
	 * @param solution
	 * @return
	 */

	private boolean canMap(Node currentNode, int index, List<Integer> edges_i,
			HashMap<Integer, Node> solution) {
		Iterator<Integer> edgeIterator = edges_i.iterator();

		// for the edges that should exist in the current node
		while (edgeIterator.hasNext()) {
			int indexNeighbor = edgeIterator.next();

			if (solution.containsKey(indexNeighbor)) {
				Node node_index = solution.get(indexNeighbor);
				Iterator<Relationship> currentNeighbors = currentNode
						.getRelationships().iterator();
				// if the current node is connected with Node_index
				boolean edgeExists = false;
				//one of the neighbor - matched 
				while (currentNeighbors.hasNext()) {
					if (currentNeighbors.next().getEndNode().getId() == node_index.getId()) {
						// current edge is matched.
						edgeExists = true;
					}
				}
				// if one of the edge that has to be present is not present =>
				if (!edgeExists)
					return false;

			}
		}
		return true;
	}

	
	public void printSearchSpace(){
		int index=0;
		for (List<Node> nodeBucket : nodeBuckets) {
			int numberOfNodes = nodeBucket.size();
			// update the least search space
			String label = nodeBucket.get(0).getLabels().iterator().next().name();
			List<String> profile = this.nodeProfileMap.get(index);
			System.out.println("No of nodes present in Bucket :"+index+" is "+numberOfNodes+"  with label "+label+" with profile"+profile);
			index++;
		}
	
	}
	
	
	
	public static void main(String[] args) {
		if (args.length != 4) {
			System.out
					.println("Enter the path of  target database, path of the query graph file, igraph/protein as 1 or 0");
			return;
		}
		NaiveGraphMatching naiveMatcher = new NaiveGraphMatching();

			//SAMPLE - DATA
			//File  igraphOrProtein = new File("C:\\Users\\laksh\\Documents\\Graph_Databases\\Assignments\\A4\\Proteins\\Proteins\\Proteins\\query\\backbones_1EMA.8.sub.grf");
			//String proteinDb = "C:\\Users\\laksh\\Documents\\Graph_Databases\\Neo4j\\DB\\backbones_1O54.grf\\";
			String db = args[0];
			File igraphOrProtein = new File(args[1]);
			
			naiveMatcher.graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(new File(db));
			try (Transaction tx = naiveMatcher.graphDb.beginTx()) {

			System.out.println("Loading db from the path..." + db);
			System.out.println("Loading query from the path..." + igraphOrProtein.getAbsolutePath());

					
			
			/**
			 * PREPROCESSING OF QUERY - FROM THE QUERY FILE  IDENTIFY THE GRAPH AND ITS PROFILES
			 */
			// get the query from the text file mentioned
			//naiveMatcher.getNodesEdgesIgraphQuery(iGraphQuery);
			System.out.println("Selected file-type to process..."+ args[2]);
			if(args[2].equals("1")){
			naiveMatcher.getNodesEdgesProteinQuery(igraphOrProtein);
			}else{
				naiveMatcher.getNodesEdgesIgraphQuery(igraphOrProtein);

			}
			// set the edges and profiles
			// set the potential-nodes in the corresponding bucket of nodes
			ResourceIterator<Node>[] filteredNodes = naiveMatcher.getNodes();
			naiveMatcher.setIncomingEdges();
			naiveMatcher.nodeBuckets = naiveMatcher.filterProfileBasedNodes(filteredNodes);
			naiveMatcher.printQuery();
			naiveMatcher.printSearchSpace();
			// compute the processing order		
			if(args[3].equals("1")){
				System.out.println("\n\nComputing processing order to execute the graph matching");
			int[] processingOrder=naiveMatcher.computeProcessingOrder();
			
			// start the graph matching and process in the order of processing order
			naiveMatcher.startGraphMatching(0, new HashMap<Integer, Node>(), processingOrder);
			}
			else{
				naiveMatcher.startGraphMatching(0, new HashMap<Integer, Node>(), null);

			}
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	/***
	 * Debug utility to print node names and id to check 
	 */
	private void printQuery(){
		Iterator<String> nodeName =this.nodeBucketIdMap.keySet().iterator();
		//print all nodes
		while(nodeName.hasNext()){
			String nodename = nodeName.next();
			System.out.println("Name:"+nodename+"\t index: "+ nodeBucketIdMap.get(nodename)+"\t Label:"+this.nodeLabelMap.get(nodename));
		}
		//print the incoming edges to every node
		int Fromindex=0;
		for(Integer[] edges : this.incomingEdges){
				for(Integer toIndex : edges){
					System.out.println(this.incomingEdges[Fromindex][0]);
					if(incomingEdges[Fromindex][toIndex] ==1){
					System.out.println("From : "+Fromindex+"\t"+" toIndex:"+toIndex);
					}
				}
				Fromindex++;
		}
	}
	
	private int[] computeProcessingOrder() {
		// from the query node compute the processing order
		// the order is the index of a node-bucket - we have N nodes in query we
		// have
		int[] processingOrder = new int[nodeBuckets.length];
		Set<Integer> processedIndices = new LinkedHashSet();
		int nodeIndex = -1;
		long leastSearchSpace = Long.MAX_VALUE;
		// identify the least search space node
		int index = 0;
		for (List<Node> nodeBucket : nodeBuckets) {
			int numberOfNodes = nodeBucket.size();
			// update the least search space
			if (leastSearchSpace > numberOfNodes) {
				leastSearchSpace = numberOfNodes;
				nodeIndex = index;
			}
			index++;
		}
		// add the current index to the set processed indices
		// if U-0 indexed at 2, is processed,
		processingOrder[0] = nodeIndex;
		processedIndices.add(nodeIndex);
		List<Integer> indicesToProcess = new LinkedList<Integer>();
		for (int i = 0; i < nodeBuckets.length; i++) {
			if (!processedIndices.contains(i)) {
				indicesToProcess.add(i);
			}
		}

		processingOrder = computeProcessingOrderUtil(processingOrder, 1,
				processedIndices, indicesToProcess);
		return processingOrder;
	}

	/**
	 * the recursive way of finding the other processing order
	 * it returns the array [1,3,4,52] => which is the order in which the nodes have to be processed
	 * This means u[1], u[3].. is the order which graph will be matched
	 * @param existingOrder
	 * @param startIndex
	 * @return
	 */
	private int[] computeProcessingOrderUtil(int[] existingOrder,
			int currentIndex, Set processedIndices,
			List<Integer> indicesToProcess) {
		
		//if all the indices have been processed
		if(currentIndex == nodeBuckets.length)
			return existingOrder;
		
		// incoming edges to a node
		Iterator<Integer> indexIter = indicesToProcess.iterator();

		// for every node that has not been processed pick the one with the
		// least cost of joining -currently
		double leastCostSoFar = Double.MAX_VALUE;
		int leastCostingIndex =-1;
		while (indexIter.hasNext()) {

			int node_i = indexIter.next();
			Integer[]edges_i = this.incomingEdges[node_i];
			List<Integer> edges = new LinkedList<Integer>();
			for(int i=0;i<edges_i.length;i++){
				if(edges_i[i] ==1){
					edges.add(i);
				}
			}
			List<Integer> incoming_i =  edges;
			// how many elements in list incoming_i is present in
			// processedIndices is the value of edge_present
			Iterator<Integer> edgeIterator = incoming_i.iterator();
			int edge_present = 0;
			while (edgeIterator.hasNext()) {
				int start_node_index = edgeIterator.next();
				if (processedIndices.contains(start_node_index)) {
					edge_present++;
				}
			}
			int search_space_i = nodeBuckets[node_i].size();
			double currentNodeCost =search_space_i * Math.pow(0.5, edge_present);
			if(currentNodeCost < leastCostSoFar){
				leastCostSoFar= currentNodeCost;
				leastCostingIndex = node_i;
			}
		}

		//after processing all the nodes for current iteration
		existingOrder[currentIndex] = leastCostingIndex;
		indicesToProcess.remove(new Integer(leastCostingIndex));
		processedIndices.add(leastCostingIndex);

		return computeProcessingOrderUtil(existingOrder,currentIndex+1, processedIndices,indicesToProcess);
	}

}
