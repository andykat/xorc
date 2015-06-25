using System;
using System.Collections;
using System.Collections.Generic;
namespace xorc
{
	class MainClass
	{
		public int inputN = 2;
		public int outputN = 1;
		public int populationSize = 200;

		//number of species limits
		public int maxSpeciesN = 30;
		public int minSpeciesN = 5;

		//number of allowable generations for protecting innovation through speciation
		public int allSpeciesLastLimit = 40;
		public int speciesLastLimit = 30;
		public double allSpeciesImprovementThreshold = 0.001;
		public double speciesImprovementThreshold = 0.001;

		//species classification
		public double weightWeight = 0.4;
		public double excessWeight = 1.0;
		public double disjointWeight = 1.0;
		//classification differences are adjusted to keep the number of different Species stable.
		public double classificationDifferenceThreshold = 1.0;
		public double currentCDT = 1.0;
		public double cdtMultiplier = 0.1;

		//competition
		public double speciesKeep = 0.4; //percent of networks in a species that are not killed each generation. These species are also the ones that breed the new generation.

		//all different types of mutation chances
		public double mutateConnectionsChance = 0.25; //??
		public double perturbChance = 0.75;
		public double perturbPercent = 0.2;
		public double crossoverChance = 0.75;
		public double linkMutationChance = 0.05;
		public double nodeMutationChance = 0.05;
		public double disableMutationKeep = 0.75;
		public double enableMutationChance = 0.2;



		/// ///////////////////////////////////////////////////////////////////
		/// /// ///////////////////////////////////////////////////////////////
		/// /// ///////////////////////////////////////////////////////////////
		/// /// /// //////////////////Array Variable///////////////////////////
		/// /// /// //////////////////Initializations//////////////////////////
		/// /// /// ///////////////////////////////////////////////////////////
		/// /// /// ///////////////////////////////////////////////////////////
		/// /// /// ///////////////////////////////////////////////////////////
		/// /// /// ///////////////////////////////////////////////////////////


		public List<Network> networks;
		public List<int> innovationInputs;
		public List<int> innovationOutputs;
		public int neuronN;
		//the first genome in a species list is the one used to compare for new genomes
		public List<List<int>> species;



		public List<double> speciesLastMaxFitness;
		public List<int> speciesLastFitnessImprovement;

		public double allSpeciesLastMaxFitness = -9999.0;
		public int allSpeciesLastFitnessImprovement = 0;

		public int startingIndexOfNewGenomes = 0;

		public Random rand = new Random();
		public static void Main (string[] args)
		{
			//testNetwork ();
		}
		public static void populationSetup(){
			innovationInputs = new List<int> ();
			innovationOutputs = new List<int> ();
			for(int i=0;i<populationSize;i++){
				Network net = new Network();

				//add inputs
				for(int j=0;j<inputN;j++){
					Neuron n = new Neuron();
					net.addNeuron(n);
					net.addInput(j);
				}

				//add outputs
				for(int j=0;j<outputN;j++){
					Neuron n = new Neuron();
					net.addNeuron(n);
					net.addOutput(inputN + j);
				}

				//add edges to outputs
				for(int j=0;j<inputN;j++){
					for(int k=0;k<outputN;k++){
						//add a random weight between -2 and 2
						bool temp = net.addEdge(j, k+inputN, hyperbolicTangent( rand.NextDouble() * 4.0 - 2.0));
					}
				}
			}

			//initialize species
			species = new List<List<int>> ();


		}

		//One simulation loop (one generation)
		public static void loop(){
			//split genomes into species
			for(int i=startingIndexOfNewGenomes;i<networks.Count;i++){
				int speciesIndex = -1;
				for(int j=0;j<species.Count;j++){
					if(classification(i, species[j][0]) < currentCDT){
						speciesIndex = j;
						break;
					}
				}

				//did not find similar species
				if(speciesIndex == -1){
					//add to new species
					List<int> tList = new List<int>();
					tList.Add(i);
					species.Add(tList);
				}
				else{
					species[speciesIndex].Add(i);
				}
			}

			//adjust classificationDifferentThreshold
			if (species.Count < minSpeciesN) {
				currentCDT *= (1.0 + cdtMultiplier);
			}
			if (species.Count > maxSpeciesN) {
				currentCDT *= (1.0 - cdtMultiplier);
			}

			//calculate networks' fitness
			for (int i=0; i<networks.Count; i++) {
				networks[i].fitness = xorFitness(i);
			}

			//order the fitness of networks within species
			for(int i=0;i<species.Count;i++){
				for(int j=0;j<species[i].Count;j++){
					double max = -9999.0;
					int maxIndex = 0;
					for(int k = j;k<species[i].Count;k++){
						if(networks[species[i][k]].fitness > max){
							max = networks[species[i][k]].fitness;
							maxIndex = k;
						}
					}
					//swap max with bottom
					int temp = species[i][j];
					species[i][j] = species[i][maxIndex];
					species[i][maxIndex] = temp;
				}
			}

			int remainingNumberOfGenomes = 0;
			//delete poor performing networks in a species
			for (int i=0; i<species.Count; i++) {
				int speciesDeleted = (int)(((double) species[i].Count) * (1.0 - speciesKeep));
				for (int j=0; j<speciesDeleted; j++) {
					species[i].RemoveAt(species[i].Count - 1);
				}
				remainingNumberOfGenomes += species[i].Count;
			}


			double populationAverageFitness = 0.0;
			//calculate the average fitness of species
			List<int> speciesAverageFitness = new List<int> ();
			for (int i=0; i<species.Count; i++) {
				double sum = 0.0;
				for(int j=0;j<species[i].Count;j++){
					sum += networks[species[i][j]].fitness;
				}
				populationAverageFitness += sum;
				speciesAverageFitness.Add(sum / ((double)species[i].Count));
			}

			populationAverageFitness /= ((double)networks.Count);

			if ((populationAverageFitness / allSpeciesLastMaxFitness) > 1.0 + allSpeciesImprovementThreshold) {
				allSpeciesLastFitnessImprovement = 0;
				allSpeciesLastMaxFitness = populationAverageFitness;
			} else {
				allSpeciesLastFitnessImprovement++;
			}

			//average population has not improved for a while, kill off all but top 2 species
			if (allSpeciesLastFitnessImprovement > allSpeciesLastLimit) {
				for(int i=0;i<species.Count - 2;i++){
					species.RemoveAt(species.Count - 1);
				}
				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				//remove all the species last fitness
				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			}

			//get ranking of species
			int speciesRanking = new int[species.Count];
			for (int i=0; i<species.Count; i++) {
				speciesRanking [i] = -1;
			}
			for (int i=0; i<species.Count; i++) {
				double maxSpeciesFitness = -9999.0;
				double maxSpeciesIndex = -1;
				for (int j=0; j<species.Count; j++) {
					if(speciesRanking[j] == -1){
						if(speciesAverageFitness[j] > maxSpeciesFitness){
							maxSpeciesFitness = speciesAverageFitness[j];
							maxSpeciesIndex = j;
						}
					}
				}
				speciesRanking[maxSpeciesIndex] = species.Count - i;
			}

			double numberOfChildren = ((double)populationSize - ((double) remainingNumberOfGenomes));
			List<Network> childrenGenomes = new List<Network> ();
			//breed new genomes
			for (int i=0; i<species.Count; i++) {
				int speciesChildrenN = getChildrenN((double)speciesRanking[i], (double)species.Count, numberOfChildren);
				for(int j=0;j<speciesChildrenN;j++){
					int parent0 = getParent(species[i].Count);
					int parent1 = getParent(species[i].Count); 

					List<Edge> childEdges;
					if(parent0 < parent1){
						childEdges = breed(species[i][parent0], species[i][parent1]);
					}
					else{
						childEdges = breed(species[i][parent1], species[i][parent0]);
					}


					//add new genome to childrenGenomes
				}
			}


			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			//remove deleted genomes and update species list
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			List<Network> parentGenomes = new List<Network> ();
			List<List<int>> newSpecies = new List<List<int>> ();
			for (int i=0; i<species.Count; i++) {
				List<int> tList = new List<int>();
				newSpecies.Add(tList);
				for(int j=0;j<species[i].Count;j++){
					parentGenomes.Add(networks[species[i][j]]);
					newSpecies[i].Add(parentGenomes.Count-1);
				}
			}
			species = newSpecies;
			networks = parentGenomes;
			startingIndexOfNewGenomes = networks.Count;
			for (int i=0; i<childrenGenomes; i++) {
				networks.Add (childrenGenomes [i]);
			}
			
		}

		private List<Edge> breed(int parentA, int parentB)
		{
			List<Edge> childEdges = new List<Edge>();
			int AIndex = 0;
			int BIndex = 0;

			while (AIndex < networks[parentA].edges.Count && BIndex < networks[parentB].edges.Count) {
				if(networks[parentA].edges[AIndex].innovation == networks[parentB].edges[BIndex].innovation){
					if(rand.Next(0,2) == 0){
						childEdges.Add (networks[parentA].edges[AIndex]);
					}
					else{
						childEdges.Add (networks[parentB].edges[BIndex]);
					}
					AIndex++;
					BIndex++;
				}
				else if(networks[parentA].edges[AIndex].innovation > networks[parentB].edges[BIndex].innovation){
					BIndex++;
				}
				else{
					childEdges.Add (networks[parentA].edges[AIndex]);
					AIndex++;
				}
			}
			//excess edges left over from dominant parent
			if (BIndex == networks [parentB].edges.Count) {
				for (int i=AIndex; i<networks[parentA].edges.Count; i++) {
					childEdges.Add (networks [parentA].edges [i]);
				}
			}

			//perturb it
			for (int i=0; i<childEdges.Count; i++) {
				double rand01 = rand.NextDouble();
				//perturb it! sore aru
				if(rand01 < perturbChance){
					childEdges[i].weight = childEdges[i].weight * (1.0 + (rand.NextDouble()*perturbPercent * 2.0 - perturbPercent));
				}
			}
			
			//new links
			double rand2 = rand.NextDouble();
			if (rand2 < linkMutationChance) {

			}
			
			//new node
			rand2 = rand.NextDouble ();
			if (rand2 < nodeMutationChance) {

			}

			return childEdges;
		}


		//returns how far apart the two genomes are
		private int classification(int netIndexA, int netIndexB){
			int AIndex = 0;
			int BIndex = 0;
			double disjointN = 0.0;
			double excessN = 0.0;
			double weightDifference = 0.0;
			while (AIndex < networks[netIndexA].edges.Count && BIndex < networks[netIndexB].edges.Count) {
				if(networks[netIndexA].edges[AIndex].innovation == networks[netIndexB].edges[BIndex].innovation){
					weightDifference += Math.Abs(networks[netIndexA].edges[AIndex].weight - networks[netIndexB].edges[BIndex].weight);
					AIndex++;
					BIndex++;
				}
				else if(networks[netIndexA].edges[AIndex].innovation > networks[netIndexB].edges[BIndex].innovation){
					disjointN += 1.0;
					BIndex++;
				}
				else{
					disjointN += 1.0;
					AIndex++;
				}
			}
			if (AIndex == networks [netIndexA].edges.Count) {
				excessN = networks [netIndexB].edges.Count - BIndex;
			} else if (BIndex == networks [netIndexB].edges.Count) {
				excessN = networks [netIndexA].edges.Count - AIndex;
			} else {
				Console.WriteLine("network classification error");
			}

			return ((weightWeight * weightDifference) + (disjointWeight * disjointN + excessWeight * excessN) / bigger (networks [netIndexA].edges.Count, networks [netIndexB].edges.Count));
		}

		//keeps indexes on the innovation count of each gene. New genes increase innovation
		private int getInnovation(int input, int output){
			for (int i=0; i<innovationInputs.Count; i++) {
				if(input == innovationInputs[i]){
					if(output == innovationOutputs[i])
					{
						return i;
					}
				}
			}
			innovationInputs.Add (input);
			innovationOutputs.Add (output);
			return (innovationInputs.Count - 1);
		}

		//calculates the different between the ideal answer from all four possible cases
		private double xorFitness(int index)
		{
			List<double> t0 = new List<double> ();
			t0.Add (0.0);
			t0.Add (0.0);
			List<double> t1 = new List<double> ();
			t1.Add (1.0);
			t1.Add (0.0);
			List<double> t2 = new List<double> ();
			t2.Add (0.0);
			t2.Add (1.0);
			List<double> t3 = new List<double> ();
			t3.Add (1.0);
			t3.Add (1.0);
			double fitness = xorFitnessSingle(0.0, 0.0, networks[index].calculateOutput(t0)) + 
							 xorFitnessSingle(1.0, 0.0, networks[index].calculateOutput(t1)) +
							 xorFitnessSingle(0.0, 1.0, networks[index].calculateOutput(t2)) + 
							 xorFitnessSingle(1.0, 1.0, networks[index].calculateOutput(t3));
			return fitness;
		}

		private int getParent(int total){
			int totalRandoms = 0;
			if (total % 2 == 0) {
				totalRandoms = (total / 2) * (total + 1);
			} else {
				totalRandoms = ((total+1)/2) * total;
			}
			totalRandoms += totalRandoms / 4;
			int randomInt = rand.Next (1, totalRandoms+1);
			int sum = 0;
			for (int i=1; i<total+1; i++) {
				sum += i;
				if (randomInt <= i) {
					return (total - i);
				}
			}
			return 0;
		}

		private int getChildrenN(double ranking, double total, double totalChildren){
			return (int)(ranking / (total * (total + 1.0) / 2.0) * totalChildren);
		}


		//calculates the fitness of a single xor case
		private double xorFitnessSingle(double inputA, double inputB, double output)
		{
			double answer = inputA + inputB;
			if (answer > 1.01) {
				answer = 0.0;
			}
			return 1.0 / (1.0 + Math.Abs (answer - output));
		}

		private void hyperbolicTangent(double x)
		{
			double e2z = Math.Pow (Math.E, 2.0 * x);
			return (e2z - 1.0) / (e2z + 1.0);
		}

		private int bigger(int a, int b){
			if (a > b) {
				return a;
			}
			return b;
		}
		private int smaller(int a, int b){
			if (a < b) {
				return a;
			}
			return b;
		}

		//deprecated
		public static void testNetwork(){
			Network a = new Network ();
			for (int i=0; i<7; i++) {
				Neuron n = new Neuron ();
				a.addNeuron (n);
			}
			
			a.addInput (0);
			a.addInput (1);
			a.addInput (2);
			a.addOutput (3);
			bool x = a.addEdge (1, 4, 0.1);
			x = a.addEdge (2, 4, 0.2);
			x = a.addEdge (1, 5, 0.3);
			x = a.addEdge (0, 5, 0.4);
			x = a.addEdge (4, 5, 0.5);
			x = a.addEdge (5, 6, 0.6);
			x = a.addEdge (6, 3, 0.7);
			
			List<double> inputs = new List<double> ();
			inputs.Add (-1);
			inputs.Add (0.5);
			inputs.Add (1.0);
			List<double> outputs = a.calculateOutput (inputs);
			for (int i=0; i<outputs.Count; i++) {
				Console.WriteLine ("i:" + outputs [i]);
			}
			Console.WriteLine ("done");
		}
	}
}
