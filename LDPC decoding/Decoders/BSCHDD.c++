#include <bits/stdc++.h>
#include <ctime>
using namespace std;
// Declaration of Node structure
// Every Node in the Tanner graph is represented by this structure.
// It has a value which stores the current value the VN or CN has.
// It has a array Message_Recived which stores the messages recived by the nodes connected to it
// An array of pointers which will point to the nodes where they have a connnection
// If VN1 has a connection with cn1,cn2, cn3, then the pointers in this connected_nodes array will point to
// a position in Message_recived array of cn1,cn2,cn3
// When sending messages, it will directly accesss the position and changes the value which is to be sent.
// In this way, the messages are passed between nodes
struct Node
{
  int value;
  int **connected_nodes;
  int *Message_Recieved;
};
// A Function for introducing  the errors which will act like a Binary Symmetric channel
// It will take a transmitted array, its size, and a float value indicating error probability
// It will create a new array by adding errors with given probability and sends the message containing error
int *Bsc_channel(int *transmitted, int n, float p)
{
  int *Recieved = new int[n];
  for (int i = 0; i < n; i++)
  {

    if (rand() % 100 < p * 100)
      Recieved[i] = !transmitted[i];
    else
      Recieved[i] = transmitted[i];
  }

  return Recieved;
}
// A boolean function which takes the array of nodes, especially check Nodes and its size.
// As each node has a value, check nodes also has a value indicating the parity of the nodes connected to it
// So, to satisfy parity, they should be all zeros
// this function will return true if all the parities for each check node are matched
// it means all the check nodes should have a value 0, else it will return false
bool are_Satisfied(Node *cns, int n)
{
  for (int i = 0; i < n; i++)
  {
    if (cns[i].value != 0)
      return false;
  }
  return true;
}
// A boolean function are_matched
// It takes the array of vns as they contain the decode bit and the transmited message and compares
// each value, and also the size of the message is passed to this function
// for every operation we are passing the size so that the array indexes should not be out of bound
bool are_matched(Node *vns, int *Transmitted, int n)
{
  for (int i = 0; i < n; i++)
  {
    if (vns[i].value != Transmitted[i])
      return false;
  }
  return true;
}
// A function named majority, which will do majority decoding
// It takes an array of messages, and a value which is recived from the channel and an integer for
// size of the array,
// specifically, the array contains messages from cn to vn and it will do majority decoding
// between the ri(recieved from channel ) and these messages
// It will return 1 if more 1's are there than 0's and return 0's if more zeros are there
int majority(int *arr, int v, int n)
{
  int n1 = 0, n0 = 0;
  for (int i = 0; i < n; i++)
  {
    if (arr[i] == 1)
      n1++;
    else
      n0++;
  }
  v == 1 ? n1++ : n0++;
  return n1 > n0 ? 1 : 0;
}
int main(int argc, char const *argv[])
{
  ifstream File;
  File.open("HMATRIX_3792x5056.txt");
  // Now with this matrix, we can calculate the number of Check nodes and variable nodes
  int Number_of_CNs;
  int Number_of_VNs;
  File >> Number_of_CNs >> Number_of_VNs;
  cout << Number_of_CNs << " " << Number_of_VNs << endl;
  int **H = new int *[Number_of_CNs];
  for (int i = 0; i < Number_of_CNs; i++)
    H[i] = new int[Number_of_VNs];
  for (int i = 0; i < Number_of_CNs; i++)
    for (int j = 0; j < Number_of_VNs; j++)
    {
      if (!File.eof())
        File >> H[i][j];
    }
  File.close();
  int dc[Number_of_CNs], dv[Number_of_VNs];
  // two arrays dc,dv to mantain degree for each check node and variable node
  // Next we will go trough each row and calculate the degree of each check nodes
  for (int i = 0; i < Number_of_CNs; i++)
  {
    int sum = 0;
    for (int j = 0; j < Number_of_VNs; j++)
      sum += H[i][j];
    dc[i] = sum;
  }
  // similarly for each column calcuating the degree of each VN
  for (int i = 0; i < Number_of_VNs; i++)
  {
    int sum = 0;
    for (int j = 0; j < Number_of_CNs; j++)
      sum += H[j][i];
    dv[i] = sum;
  }
  // Now creating an array of nodes
  // ONe for Checkn Nodes and one for variable noddes
  Node *Check_Nodes = new Node[Number_of_CNs];
  Node *Variable_Nodes = new Node[Number_of_VNs];
  // Next intialising the arrays in each node with the size of its degree
  for (int i = 0; i < Number_of_CNs; i++)
  {
    Check_Nodes[i].connected_nodes = new int *[dc[i]];
    Check_Nodes[i].Message_Recieved = new int[dc[i]];
    for (int j = 0; j < dc[i]; j++)
      Check_Nodes[i].Message_Recieved[j] = 4;
    // Intial value 4 is given indicating that the connection are not made for this
  }
  for (int i = 0; i < Number_of_VNs; i++)
  {
    Variable_Nodes[i].connected_nodes = new int *[dv[i]];
    Variable_Nodes[i].Message_Recieved = new int[dv[i]];
    for (int j = 0; j < dv[i]; j++)
      Variable_Nodes[i].Message_Recieved[j] = 4;
    // 4 indicating it is not connected to any other node
  }
  // Now making the connections
  for (int i = 0; i < Number_of_CNs; i++)
  {
    // For each check node going through all the vns and
    // if there is 1 in the H matrix in these position then there is a connnection to it
    // Next , checking the position in Message_recieved array for that vn where ther is no connnection
    // It means it should be 4 as we intialised with 4
    int j = 0;
    for (int k = 0; k < Number_of_VNs; k++)
    {
      if (H[i][k] != 1)
        continue;
      Node *vn = Variable_Nodes + k;
      int l;
      for (l = 0; l < dv[k]; l++)
      {
        if (vn->Message_Recieved[l] == 4) // if 4 is found then we can make connection
          break;
      }
      // connect both the cn and vn and change the value in that position from 4 to 3 indicating that
      // the connnection is made
      Check_Nodes[i].connected_nodes[j] = &vn->Message_Recieved[l];
      *Check_Nodes[i].connected_nodes[j] = 3;
      j++;
      vn = NULL;
      delete vn;
    }
  }
  for (int i = 0; i < Number_of_VNs; i++)
  {
    // similarly as in cns,make connections for vns
    // search the position to be connected and change it to 3
    int j = 0;
    for (int k = 0; k < Number_of_CNs; k++)
    {
      if (H[k][i] != 1)
        continue;
      Node *cn = Check_Nodes + k;
      int l;
      for (l = 0; l < dc[k]; l++)
      {
        if (cn->Message_Recieved[l] == 4)
          break;
      }
      Variable_Nodes[i].connected_nodes[j] = &cn->Message_Recieved[l];
      *Variable_Nodes[i].connected_nodes[j] = 3;
      j++;
      cn = NULL;
      delete cn;
    }
  }
  // Until this part, we created check nodes, intialised with the degrees and
  // made the connections between check nodes and variable nodes
  // with the help of pointers so that we can directly send the message
  // From now the algorithm will start
  int Transmitted[Number_of_VNs] = {0};             // intialising all-zero codeword
  int length = 101;                                // length of the p vector we are plotting from 0 to 1
  float pVec[length] = {0};                         // p vector for error probaility
  int MAXdc = *max_element(dc, dc + Number_of_CNs); // finding maximum value of dc and dv
  int MAXdv = *max_element(dv, dv + Number_of_VNs);
  //The reason for findng Maxdc and Maxdv is
  // To allocate some arrays for handling messages
  // As we are doing through arrays and use indexes everywhere,there will be no problem if the size of array is
  // more
  // So , that we can only allocate memory one time, and use these array everywhere
  // Otherwise memory related problems will occur for large H matrices
  // and also the program becomes slow if we allocate memory each time
  int *Message_from_dcminus1nodes = new int[MAXdc];
  int *messages = new int[MAXdc > MAXdv ? MAXdc : MAXdv];
  int *messageFromDvminus1Nodes = new int[MAXdv];
  Node *vn = NULL, *cn = NULL; // creating a two node pointers, to use temporarily
  float pSuccessVec[length] = {0};
  // A vector for storing the succesfull decoding probabilty
  // Now creating a pvector from 0 to 1
  float intial_value = 0, increment = 1 / (float(length - 1));
  for (int i = 0; i < length; i++)
  {
    pVec[i] = intial_value;
    intial_value += increment;
  }
  // Then for each value of p run Nsim times
  for (int pe = 0; pe < length; pe++)
  {
    int Nsim = 1000; // Number of simulations
    int Nerr = 0;
    int Ncorr = 0;
    for (int num = 1; num <= Nsim; num++)
    {
      int *Recieved = Bsc_channel(Transmitted, Number_of_VNs, pVec[pe]);
      // Passing the message into BSc channel and storing it in recieved array
      for (int i = 0; i < Number_of_VNs; i++)
      {
        // In the intial iteration ,update the values of vns with the recieved message and send
        // these values directly to the cns
        Node *vn = Variable_Nodes + i;
        vn->value = Recieved[i];
        for (int j = 0; j < dv[i]; j++)
          *vn->connected_nodes[j] = vn->value;
        // by using the address of the location we can directly send messages
      }
      int iteration_index = 0; // Intialisng iteration index to 0
      int max_iterations = 50;  // maximum allowed iterations
      while (iteration_index <= max_iterations)
      {
        iteration_index++;
        // CN to VN message passing
        for (int i = 0; i < Number_of_CNs; i++)
        {
          cn = Check_Nodes + i;
          // collecting messages recieved in an array
          for (int j = 0; j < dc[i]; j++)
            messages[j] = cn->Message_Recieved[j];
          // sending messages to all Vns connected to it
          for (int j = 0; j < dc[i]; j++)
          {
            int sum = 0;
            for (int k = 0; k < dc[i]; k++)
            {
              if (k != j) // sending messages to all Vns connected to it
                sum += messages[k];
            }
            *cn->connected_nodes[j] = sum % 2;
            // sending the mod 2 sum to vns
          }
          cn->value = accumulate(messages, messages + dc[i], 0) % 2;
          // updating the cn value as modulo two sum of all messages connnected to it
        }
        // If all the check nodes are satisfied then terminate
        if (are_Satisfied(Check_Nodes, Number_of_CNs))
          break;
        // VN to CN message passing
        for (int i = 0; i < Number_of_VNs; i++)
        {
          vn = Variable_Nodes + i;
          for (int j = 0; j < dv[i]; j++)
            messages[j] = vn->Message_Recieved[j];
          // collecting messages
          for (int j = 0; j < dv[i]; j++)
          {
            int counter = 0;
            for (int k = 0; k < dv[i]; k++)
            {
              if (k != j)
                messageFromDvminus1Nodes[counter++] = messages[k];
            }
            *vn->connected_nodes[j] = majority(messageFromDvminus1Nodes, Recieved[i], dv[i] - 1);
            // sending the majority voting value to each cn
          }
          vn->value = majority(messages, Recieved[i], dv[i]);
          // Updating the value of vn
        }
      }
      // After decoding, if it is matchin with transmitted then increase Ncorr
      if (are_matched(Variable_Nodes, Transmitted, Number_of_VNs))
        Ncorr++;
      else
        Nerr++;
    }
    cout<<pe;
    pSuccessVec[pe] = (float(Ncorr) / Nsim);
    // storing the value in array for succesfull probability
  }

  // Plotting graph
  FILE *fp = NULL;
  FILE *gnupipe = NULL;
  char *GnuCommands[] = {"set title \"Successful Decoding for BSC Hard decision\"",
                         "set key noautotitle",
                         "set xlabel \"BSC error probability(p)\"",
                         "set ylabel \"Successful decoding probability\"",
                          "plot 'data.tmp' with lines linecolor 2 linewidth 2",
                          "set terminal wxt size 700,600",
                          "set grid",
                          "replot"};
  fp = fopen("data.tmp", "w");
  gnupipe = _popen("gnuplot -persitent,", "w");
  for (int i = 0; i < length; i++)
  {
    fprintf(fp, "%f %f\n", pVec[i], pSuccessVec[i]);
  }
  for (int i = 0; i < 8; i++)
  {
    fprintf(gnupipe, "%s\n", GnuCommands[i]);
  }
  cout << "Done,By this time the graph should be visible to you";
  return 0;
}
