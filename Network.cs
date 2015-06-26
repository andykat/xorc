//------------------------------------------------------------------------------
// <auto-generated>
//     This code was generated by a tool.
//     Runtime Version:4.0.30319.34209
//
//     Changes to this file may cause incorrect behavior and will be lost if
//     the code is regenerated.
// </auto-generated>
//------------------------------------------------------------------------------
using System;
using System.Collections;
using System.Collections.Generic;
namespace xorc
{
	public class Network{
		public List<Neuron> neurons;
		public List<int> inputNeuronIndexes;
		public List<int> outputNeuronIndexes;
		public List<Edge> edges;
		public double fitness = 0.0;
		private int finishedNeuronCounter = -9999;
		public Network(){
			neurons = new List<Neuron> ();
			inputNeuronIndexes = new List<int>();
			outputNeuronIndexes = new List<int> ();
			edges = new List<Edge>();
		}
		
		public void addNeuron(Neuron n){
			neurons.Add (n);
		}
		public void addInput(int index){
			inputNeuronIndexes.Add (index);
		}
		public void addOutput(int index){
			outputNeuronIndexes.Add (index);
		}
		
		public bool addEdge(int input, int output, double weight, int innovation){
			//check if edge is already added
			for (int i=0; i<neurons[input].outputEdges.Count; i++) {
				if (output == neurons [input].outputEdges [i].outNeuron) {
					return false;
				}
			}

			//find if this edge leads to a cycle
			List<int> queue = new List<int> ();
			queue.Add (output);
			while (queue.Count>0) {
				int cIndex = queue[0];
				queue.RemoveAt(0);

				if(cIndex == input)
				{
					//creates cycle
					return false;
				}

				for(int i=0;i<neurons[cIndex].outputEdges.Count;i++){
					queue.Add(neurons[cIndex].outputEdges[i].outNeuron);
				}
			}
			
			//an input cannot be an output
			for (int i=0; i<inputNeuronIndexes.Count; i++) {
				if(output == inputNeuronIndexes[i]){
					return false;
				}
			}
			//an output cannot be an input
			for (int i=0; i<outputNeuronIndexes.Count; i++) {
				if(input == outputNeuronIndexes[i]){
					return false;
				}
			}



			//add edge
			Edge e = new Edge (input, output, weight, true, innovation);
			edges.Add (e);
			neurons [input].addOutput (e);
			neurons [output].addInput (e);

			return true;
		}

		//checks if adding/enabling this edge creates a cycle
		public bool checkEdge(int input, int output)
		{
			//find if this edge leads to a cycle
			List<int> queue = new List<int> ();
			queue.Add (output);
			while (queue.Count>0) {
				int cIndex = queue[0];
				queue.RemoveAt(0);
				
				if(cIndex == input)
				{
					//creates cycle
					return false;
				}
				
				for(int i=0;i<neurons[cIndex].outputEdges.Count;i++){
					queue.Add(neurons[cIndex].outputEdges[i].outNeuron);
				}
			}
			
			//an input cannot be an output
			for (int i=0; i<inputNeuronIndexes.Count; i++) {
				if(output == inputNeuronIndexes[i]){
					return false;
				}
			}
			//an output cannot be an input
			for (int i=0; i<outputNeuronIndexes.Count; i++) {
				if(input == outputNeuronIndexes[i]){
					return false;
				}
			}

			return true;
		}

		//calculates the output values of the network from a list of inputs
		public List<double> calculateOutput(List<double> inputs)
		{
			List<double> outputs = new List<double>();
			
			for (int i=0; i<neurons.Count; i++) {
				neurons [i].reset ();
			}
			
			//initialize inputs
			if(inputs.Count != inputNeuronIndexes.Count){
				return outputs;
			}
			//add inputs to queue
			List<int> queue = new List<int> ();
			for (int i=0; i<inputs.Count; i++) {
				neurons[inputNeuronIndexes[i]].value = inputs[i];
				queue.Add(inputNeuronIndexes[i]);
			}
			
			while (queue.Count > 0) {
				//pop
				int pIndex = queue[0];
				queue.RemoveAt(0);
				
				//add value along the output edges
				for(int i=0;i<neurons[pIndex].outputEdges.Count;i++){
					if(!neurons[pIndex].outputEdges[i].isEnabled){
						continue;
					}
					int outIndex = neurons[pIndex].outputEdges[i].outNeuron;
					neurons[outIndex].addValue(neurons[pIndex].value, neurons[pIndex].outputEdges[i].weight);
					
					//all input signals have been passed to outNeuron
					if(neurons[outIndex].inputCounter < 1){
						if(neurons[outIndex].inputCounter == 0){
							//compute value
							neurons[outIndex].computeValue();
							
							//reduce counter
							neurons[outIndex].inputCounter = finishedNeuronCounter;
							
							//add to queue
							queue.Add (outIndex);
						}
						else if(neurons[outIndex].inputCounter < finishedNeuronCounter + 1){
							Console.WriteLine("wtf added input to finished counter");
						}
						else{
							Console.WriteLine("wtf neuron with inputCounter=0 was not removed");
						}
						
					}
				}
				
			}
			
			for (int i=0; i<outputNeuronIndexes.Count; i++) {
				//output node was calculated
				if(neurons[outputNeuronIndexes[i]].inputCounter == finishedNeuronCounter){
					outputs.Add(neurons[outputNeuronIndexes[i]].value);
				}
				else
				{
					Console.WriteLine("output node was not calculated correctly");
					return outputs;
				}
			}
			for (int i=0; i<neurons.Count; i++) {
				Console.WriteLine (i + ":" + neurons [i].value);
			}
			return outputs;
		}
	}
}

