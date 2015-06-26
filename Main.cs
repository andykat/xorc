using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
namespace xorc
{
	class MainClass
	{
		public static int inputN = 2;
		public static int outputN = 1;
		public static int populationSize = 200;

		//number of species limits
		public static int maxSpeciesN = 30;
		public static int minSpeciesN = 5;

		//number of allowable generations for protecting innovation through speciation
		public static int allSpeciesLastLimit = 40;
		public static int speciesLastLimit = 30;
		public static double allSpeciesImprovementThreshold = 0.001;
		public static double speciesImprovementThreshold = 0.001;

		//species classification
		public static double weightWeight = 0.4;
		public static double excessWeight = 1.0;
		public static double disjointWeight = 1.0;
		//classification differences are adjusted to keep the number of different Species stable.
		public static double classificationDifferenceThreshold = 1.0;
		public static double currentCDT = 1.0;
		public static double cdtMultiplier = 0.1;

		//competition
		public static double speciesKeep = 0.4; //percent of networks in a species that are not killed each generation. These species are also the ones that breed the new generation.

		//all different types of mutation chances
		public static double mutateConnectionsChance = 0.25; //??
		public static double perturbChance = 0.75;
		public static double perturbPercent = 0.2;
		public static double crossoverChance = 0.75;
		public static double linkMutationChance = 0.05;
		public static double nodeMutationChance = 0.05;
		public static double disableMutationKeep = 0.75;
		public static double enableMutationChance = 0.2;



		/// ///////////////////////////////////////////////////////////////////
		/// /// ///////////////////////////////////////////////////////////////
		/// /// ///////////////////////////////////////////////////////////////
		/// /// /// //////////////////Array Variable///////////////////////////
		/// /// /// //////////////////Initializations//////////////////////////
		/// /// /// ///////////////////////////////////////////////////////////
		/// /// /// ///////////////////////////////////////////////////////////
		/// /// /// ///////////////////////////////////////////////////////////
		/// /// /// ///////////////////////////////////////////////////////////


		public static List<Network> networks;
		public static List<int> innovationInputs;
		public static List<int> innovationOutputs;
		public static int neuronN=0;
		//the first genome in a species list is the one used to compare for new genomes
		public static List<List<int>> species;



		public static List<double> speciesLastMaxFitness;
		public static List<int> speciesLastFitnessImprovement;

		public static double allSpeciesLastMaxFitness = -9999.0;
		public static int allSpeciesLastFitnessImprovement = 0;

		public static int startingIndexOfNewGenomes = 0;

		public static Random rand = new Random();
		public static void Main (string[] args)
		{
			FileStream ostrm;
			StreamWriter writer;
			TextWriter oldOut = Console.Out;
			try
			{
				ostrm = new FileStream ("./Redirect.txt", FileMode.OpenOrCreate, FileAccess.Write);
				writer = new StreamWriter (ostrm);
			}
			catch (Exception e)
			{
				Console.WriteLine ("Cannot open Redirect.txt for writing");
				Console.WriteLine (e.Message);
				return;
			}
			Console.SetOut (writer);
			Console.WriteLine ("This is a line of text");
			Console.WriteLine ("Everything written to Console.Write() or");
			Console.WriteLine ("Console.WriteLine() will be written to a file");
			Console.SetOut (oldOut);
			writer.Close();
			ostrm.Close();
			Console.WriteLine ("Done");

			//testNetwork ();
			populationSetup ();
			for (int i=0; i<100; i++) {
				Console.WriteLine("i:" + i);
				loop ();
			}

		}
		public static void populationSetup(){
			networks = new List<Network> ();
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
				neuronN = inputN + outputN;
				//add edges to outputs
				for(int j=0;j<inputN;j++){
					for(int k=0;k<outputN;k++){
						//add a random weight between -2 and 2
						bool temp = net.addEdge(j, k+inputN, hyperbolicTangent( rand.NextDouble() * 4.0 - 2.0), getInnovation(j, k+inputN));
					}
				}

				networks.Add(net);
			}

			startingIndexOfNewGenomes = 0;

			//initialize species
			species = new List<List<int>> ();
			speciesLastMaxFitness = new List<double> ();
			speciesLastFitnessImprovement = new List<int> ();
		}

		//One simulation loop (one generation)
		public static void loop(){
			//split genomes into species
			//Console.WriteLine ("networksN:" + networks.Count);
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
					speciesLastMaxFitness.Add(-9999.0);
					speciesLastFitnessImprovement.Add(0);
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
				//Console.WriteLine("nf:" + networks[i].fitness);
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
			List<double> speciesAverageFitness = new List<double> ();
			for (int i=0; i<species.Count; i++) {
				double sum = 0.0;
				for(int j=0;j<species[i].Count;j++){
					sum += networks[species[i][j]].fitness;
					//Console.WriteLine("nf:" + networks[species[i][j]].fitness);
				}
				populationAverageFitness += sum;
				//Console.WriteLine("speicescount:" + species[i].Count);
				speciesAverageFitness.Add(sum / ((double)species[i].Count));
			}

			populationAverageFitness /= ((double)networks.Count);
			Console.WriteLine ("fitness:" + populationAverageFitness);
			if ((populationAverageFitness / allSpeciesLastMaxFitness) > 1.0 + allSpeciesImprovementThreshold) {
				allSpeciesLastFitnessImprovement = 0;
				allSpeciesLastMaxFitness = populationAverageFitness;
			} else {
				allSpeciesLastFitnessImprovement++;
			}

			//average population has not improved for a while, kill off all but top 2 species
			if (allSpeciesLastFitnessImprovement > allSpeciesLastLimit) {
				for(int i=2;i<species.Count;i++){
					species.RemoveAt(i);
					speciesAverageFitness.RemoveAt(i);
					speciesLastMaxFitness.RemoveAt(i);
					speciesLastFitnessImprovement.RemoveAt(i);
					i--;
				}
			}

			for (int i=0; i<species.Count; i++) {

				//speciesLastFitnessImprovement.Add(0);
				if(speciesAverageFitness[i] > speciesLastMaxFitness[i]){
					speciesLastMaxFitness[i] = speciesAverageFitness[i];
				}
				else{
					speciesLastFitnessImprovement[i]++;
					//remove species if they fail too much
					if(speciesLastFitnessImprovement[i] > speciesLastLimit){
						species.RemoveAt(i);
						speciesAverageFitness.RemoveAt(i);
						speciesLastMaxFitness.RemoveAt(i);
						speciesLastFitnessImprovement.RemoveAt(i);
						i--;
					}
				}
			}

			//get ranking of species
			int[] speciesRanking = new int[species.Count];
			for (int i=0; i<species.Count; i++) {
				speciesRanking [i] = -1;
				//Console.WriteLine("fitnessI:" + speciesAverageFitness[i]);
			}
			//Console.WriteLine ("sc:" + species.Count + "  safc:" + speciesAverageFitness.Count);
			for (int i=0; i<species.Count; i++) {
				double maxSpeciesFitness = -9999.0;
				int maxSpeciesIndex = -1;
				for (int j=0; j<species.Count; j++) {
					if(speciesRanking[j] == -1){
						if(speciesAverageFitness[j] > maxSpeciesFitness){
							maxSpeciesFitness = speciesAverageFitness[j];
							maxSpeciesIndex = j;
						}
					}
				}
				//Console.WriteLine("msi:" + maxSpeciesIndex);
				speciesRanking[maxSpeciesIndex] = species.Count - i;
			}
			double maxFitness = 0.0;
			for (int i=0; i<species.Count; i++) {
				if (networks [species [i] [0]].fitness > maxFitness) {
					maxFitness = networks [species [i] [0]].fitness;
				}
			}
			Console.WriteLine ("maxfitness:" + maxFitness);
			double numberOfChildren = ((double)populationSize - ((double) remainingNumberOfGenomes));
			List<Network> childrenGenomes = new List<Network> ();
			//breed new genomes
			for (int i=0; i<species.Count; i++) {
				int speciesChildrenN = getChildrenN((double)speciesRanking[i], (double)species.Count, numberOfChildren);
				for(int j=0;j<speciesChildrenN;j++){
					int parent0 = getParent(species[i].Count);
					int parent1 = getParent(species[i].Count); 

					Network childNet;
					if(parent0 < parent1){
						childNet = breed(species[i][parent0], species[i][parent1]);
					}
					else{
						childNet = breed(species[i][parent1], species[i][parent0]);
					}


					//add new genome to childrenGenomes
					childrenGenomes.Add(childNet);
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
			for (int i=0; i<childrenGenomes.Count; i++) {
				networks.Add (childrenGenomes [i]);
			}
			
		}

		private static Network breed(int parentA, int parentB)
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
					childEdges[i].weight = childEdges[i].weight  + (rand.NextDouble()*perturbPercent * 2.0 - perturbPercent);
					//childEdges[i].weight = childEdges[i].weight * (1.0 + (rand.NextDouble()*perturbPercent * 2.0 - perturbPercent));
				}
			}

			//create the child network
			Network net = new Network ();
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
			for(int i=0;i<childEdges.Count;i++){
				bool flag = net.addEdge(childEdges[i].inNeuron, childEdges[i].outNeuron, childEdges[i].weight, childEdges[i].innovation);
				if(flag==false){
					Console.WriteLine("wtf add edge failed");
				}
			}


			//new links
			double rand2 = rand.NextDouble();
			if (rand2 < linkMutationChance) {
				int tCount = 0;
				while(tCount <10){
					int inputNode = rand.Next(0,neuronN - outputN);
					if(inputNode >= inputN){
						inputNode += outputN;
					}
					int outputNode = rand.Next (inputN, neuronN);
					if(net.checkEdge(inputNode, outputNode)){
						net.addEdge(inputNode, outputNode, hyperbolicTangent( rand.NextDouble() * 4.0 - 2.0), getInnovation(inputNode, outputNode));
						tCount = 20;
					}
					tCount++;
				}
			}
			
			//new node
			rand2 = rand.NextDouble ();
			if (rand2 < nodeMutationChance) {
				int tCount = 0;
				while(tCount < 10){
					int randomEdge = rand.Next(0,net.edges.Count);
					bool flag = true;
					for(int i=0;i<net.neurons[net.edges[randomEdge].inNeuron].outputEdges.Count;i++){ 
						for(int j=0;j<net.neurons[net.edges[randomEdge].outNeuron].inputEdges.Count;j++){
							if(net.neurons[net.edges[randomEdge].inNeuron].outputEdges[i] == net.neurons[net.edges[randomEdge].outNeuron].inputEdges[j]){
								flag = false;
							}
						}
					}
					//create new node
					if(flag){
						net.edges[randomEdge].isEnabled = false;
						Neuron tN = new Neuron();
						net.addNeuron(tN);
						net.addEdge(net.edges[randomEdge].inNeuron, neuronN, hyperbolicTangent( rand.NextDouble() * 4.0 - 2.0), getInnovation(net.edges[randomEdge].inNeuron, neuronN));
						net.addEdge(neuronN, net.edges[randomEdge].outNeuron, hyperbolicTangent( rand.NextDouble() * 4.0 - 2.0), getInnovation(neuronN, net.edges[randomEdge].outNeuron));
						neuronN++;
					}
					tCount++;
				}

			}

			return net;
		}


		//returns how far apart the two genomes are
		private static double classification(int netIndexA, int netIndexB){
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
		private static int getInnovation(int input, int output){
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
		private static double xorFitness(int index)
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
			//Console.WriteLine("outputs:" + networks[index].calculateOutput(t0)[0] + ""  + networks[index].calculateOutput(t1)[0] + "" + networks[index].calculateOutput(t2)[0] + "" + networks[index].calculateOutput(t3)[0] + ":"+networks[index].calculateOutput(t3).Count); 
			double fitness = xorFitnessSingle(0.0, 0.0, networks[index].calculateOutput(t0)[0]) + 
							 xorFitnessSingle(1.0, 0.0, networks[index].calculateOutput(t1)[0]) +
							 xorFitnessSingle(0.0, 1.0, networks[index].calculateOutput(t2)[0]) + 
							 xorFitnessSingle(1.0, 1.0, networks[index].calculateOutput(t3)[0]);
			return fitness;
		}

		private static int getParent(int total){
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

		private static int getChildrenN(double ranking, double total, double totalChildren){
			return (int)(ranking / (total * (total + 1.0) / 2.0) * totalChildren);
		}


		//calculates the fitness of a single xor case
		private static double xorFitnessSingle(double inputA, double inputB, double output)
		{
			double answer = inputA + inputB;
			if (answer > 1.1) {
				answer = 0.0;
			}
			return 1.0 / (1.0 + Math.Abs (answer - output));
		}

		private static double hyperbolicTangent(double x)
		{
			double e2z = Math.Pow (Math.E, 2.0 * x);
			return (e2z - 1.0) / (e2z + 1.0);
		}

		private static int bigger(int a, int b){
			if (a > b) {
				return a;
			}
			return b;
		}
		private static int smaller(int a, int b){
			if (a < b) {
				return a;
			}
			return b;
		}
		
	}
}
