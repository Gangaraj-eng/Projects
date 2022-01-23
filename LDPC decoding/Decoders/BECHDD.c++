#include <bits/stdc++.h>
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
// A Function for introducing  the erasures which will act like a Binary Erasure channel
// It will take a transmitted array, its size, and a float value indicating error probability
// It will create a new array by adding erasures with given probability and sends the message containing error
int *BEC_channel(int *transmitted, int n, float p)
{
  int *Recieved = new int[n];
  for (int i = 0; i < n; i++)
  {

    if (rand() % 100 < p * 100)
      Recieved[i] = -1;
    else
      Recieved[i] = transmitted[i];
  }
  return Recieved;
}
// A boolean function which takes the array of nodes, especially variable Nodes and its size.
// As each node has a value, variable nodes also has a value indicating the bit it has currently
// If all the bits are non-erasures it will return true
bool are_Decoded(Node *vns, int n)
{
  for (int i = 0; i < n; i++)
  {
    if (vns[i].value == -1)
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
// A function which will take an array of messages, a integer for current value of vn, and the size of array
// If all the messages are erasures, then it will return the erasure
//otherwise it will return the non_erased bit as all the non_erased bits will be equal
int Non_Erased(int *arr, int v, int n)
{
  int nnonerasures = 0, nerasures = 0, non_erased_bit;
  for (int i = 0; i < n; i++)
  {
    if (arr[i] == -1)
      nerasures++;
    else
    {
      return arr[i];
    }
  }
  if (v == -1)
    return -1;
  else
    return v;
}
// A functin pairty_check which will take an array and its size
// It will return erasure if there is any erasure present in the array
// If all are non erasures, then it will return mod 2 sum as cns are parity check nodes
int parity_check(int *arr, int n)
{
  int sum = 0;
  int erasures = 0;
  for (int i = 0; i < n; i++)
  {
    if (arr[i] == -1)
      erasures++;
    else
      sum += arr[i];
  }
  return erasures != 0 ? -1 : sum % 2;
}
int main(int argc, char const *argv[])
{
  // The declaratio of H matrix
  ifstream File;
  File.open("HMATRIX_9x12.txt");
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
  // two arrays dc,dv to mantain degree for each check node and variable nodes
  int dc[Number_of_CNs], dv[Number_of_VNs];
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
  Node *vn = NULL, *cn = NULL; // temporary pointers
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
  { // For each check node going through all the vns and
    // if there is 1 in the H matrix in these position then there is a connnection to it
    // Next , checking the position in Message_recieved array for that vn where ther is no connnection
    // It means it should be 4 as we intialised with 4
    int j = 0;
    for (int k = 0; k < Number_of_VNs; k++)
    {
      if (H[i][k] != 1)
        continue;
      vn = Variable_Nodes + k;
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
      cn = Check_Nodes + k;
      int l;
      for (l = 0; l < dc[k]; l++)
      {
        if (cn->Message_Recieved[l] == 4)
          break;
      }
      Variable_Nodes[i].connected_nodes[j] = &cn->Message_Recieved[l];
      *Variable_Nodes[i].connected_nodes[j] = 3;
      j++;
    }
  }
  // Until this part, we created check nodes, intialised with the degrees and
  // made the connections between check nodes and variable nodes
  // with the help of pointers so that we can directly send the message
  // From now the algorithm will start
  int Transmitted[Number_of_VNs] = {0}; // intialising all-zero codeword
  int length = 11;                    // length of the p vector we are plotting from 0 to 1
  int MAXdc = *max_element(dc, dc + Number_of_CNs);
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
  float pVec[length] = {0};
  float pSuccessVec[length] = {0};
  // A vector for storing the succesfull decoding probabilty
  // Now creating a pvector from 0 to 1
  float intial_value = 0, increment = 1 / (float)(length - 1);
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
      int *Recieved = BEC_channel(Transmitted, Number_of_VNs, pVec[pe]);
      // Passing the message into BSc channel and storing it in recieved array
      for (int i = 0; i < Number_of_VNs; i++)
      {
        // In the intial iteration ,update the values of vns with the recieved message and send
        // these values directly to the cns
        vn = Variable_Nodes + i;
        vn->value = Recieved[i];
        for (int j = 0; j < dv[i]; j++)
          *vn->connected_nodes[j] = vn->value;
        // by using the address of the location we can directly send messages
      }
      int iteration_index = 0; // Intialisng iteration index to 0
      int max_iterations = 10; // maximum allowed iterations
      while (iteration_index <= max_iterations)
      {
        // If all the erasures are recoverd, then break
        if (are_Decoded(Variable_Nodes, Number_of_VNs))
          break;
        iteration_index++;
        // CN to VN message passing
        for (int i = 0; i < Number_of_CNs; i++)
        {
          cn = Check_Nodes + i;
          for (int j = 0; j < dc[i]; j++) // collecting messages recieved in an array
            messages[j] = cn->Message_Recieved[j];
          // sending messages to all Vns connected to it
          for (int j = 0; j < dc[i]; j++)
          {
            int counter = 0;
            for (int k = 0; k < dc[i]; k++)
            {
              if (k != j) // sending messages to all Vns connected to it
                Message_from_dcminus1nodes[counter++] = messages[k];
            }
            *cn->connected_nodes[j] = parity_check(Message_from_dcminus1nodes, dc[i] - 1);
            // sending to vn
          }
          cn->value = parity_check(messages, dc[i]);
          // updating the cn value
        }
        for (int i = 0; i < Number_of_VNs; i++)
        {
          // similarly for each vn sending message to cns onnected to it
          vn = Variable_Nodes + i;
          for (int j = 0; j < dv[i]; j++)
            messages[j] = vn->Message_Recieved[j];
          for (int j = 0; j < dv[i]; j++)
          {
            int counter = 0;
            for (int k = 0; k < dv[i]; k++)
            {
              if (k != j)
                messageFromDvminus1Nodes[counter++] = messages[k];
            }
            *vn->connected_nodes[j] = Non_Erased(messageFromDvminus1Nodes, Recieved[i], dv[i] - 1);
          }
          vn->value = Non_Erased(messages, Recieved[i], dv[i]);
          // updating value
        }
      }
      if (are_matched(Variable_Nodes, Transmitted, Number_of_VNs))
        Ncorr++;
      else
        Nerr++;
    }
    cout << pe;
    pSuccessVec[pe] = (float(Ncorr) / Nsim); // storing the success probablity
  }
  // plotting the graph
  FILE *fp = NULL;
  FILE *gnupipe = NULL;
  char *GnuCommands[] = {"set title \"Successful Decoding for BEC Hard Decision\"",
                         "set key noautotitle",
                         "set xlabel \"BEC error probability (p)\"",
                         "set ylabel \"probability of successful decoding\"",
                         "plot 'data.tmp' with lines linecolor 3 linewidth 3",
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
  cout << "Done\n,By this time, the graph should be visible to you";
  return 0;
}
