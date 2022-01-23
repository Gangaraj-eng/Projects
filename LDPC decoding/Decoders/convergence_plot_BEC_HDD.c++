#include <bits/stdc++.h>
using namespace std;
struct Node
{
  int value;
  int **connected_nodes;
  int *Message_Recieved;
};
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
bool are_Decoded(Node *vns, int n)
{
  for (int i = 0; i < n; i++)
  {
    if (vns[i].value == -1)
      return false;
  }
  return true;
}
bool are_matched(Node *vns, int *Transmitted, int n)
{
  for (int i = 0; i < n; i++)
  {

    if (vns[i].value != Transmitted[i])
      return false;
  }
  return true;
}
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
int no_of_erasures(Node *Vns,int n)
{
  int ner=0;
  for(int i=0;i<n;i++)
  {
    if(Vns[i].value==-1)
      ner++;
  }
  return ner;
}
int main(int argc, char const *argv[])
{
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
  int dc[Number_of_CNs], dv[Number_of_VNs];
  for (int i = 0; i < Number_of_CNs; i++)
  {
    int sum = 0;
    for (int j = 0; j < Number_of_VNs; j++)
      sum += H[i][j];
    dc[i] = sum;
  }
  for (int i = 0; i < Number_of_VNs; i++)
  {
    int sum = 0;
    for (int j = 0; j < Number_of_CNs; j++)
      sum += H[j][i];
    dv[i] = sum;
  }
  Node *Check_Nodes = new Node[Number_of_CNs];
  Node *Variable_Nodes = new Node[Number_of_VNs];
  Node *vn = NULL, *cn = NULL; 
  for (int i = 0; i < Number_of_CNs; i++)
  {
    Check_Nodes[i].connected_nodes = new int *[dc[i]];
    Check_Nodes[i].Message_Recieved = new int[dc[i]];
    for (int j = 0; j < dc[i]; j++)
      Check_Nodes[i].Message_Recieved[j] = 4;
  }
  for (int i = 0; i < Number_of_VNs; i++)
  {
    Variable_Nodes[i].connected_nodes = new int *[dv[i]];
    Variable_Nodes[i].Message_Recieved = new int[dv[i]];
    for (int j = 0; j < dv[i]; j++)
      Variable_Nodes[i].Message_Recieved[j] = 4;

  }
  for (int i = 0; i < Number_of_CNs; i++)
  { 
    int j = 0;
    for (int k = 0; k < Number_of_VNs; k++)
    {
      if (H[i][k] != 1)
        continue;
      vn = Variable_Nodes + k;
      int l;
      for (l = 0; l < dv[k]; l++)
      {
        if (vn->Message_Recieved[l] == 4) 
          break;
      }

      Check_Nodes[i].connected_nodes[j] = &vn->Message_Recieved[l];
      *Check_Nodes[i].connected_nodes[j] = 3;
      j++;
    }
  }
  for (int i = 0; i < Number_of_VNs; i++)
  {
   
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
  int Transmitted[Number_of_VNs] = {0}; // intialising all-zero codeword
  int length = 31;                    // length of the p vector we are plotting from 0 to 1
  int MAXdc = *max_element(dc, dc + Number_of_CNs);
  int MAXdv = *max_element(dv, dv + Number_of_VNs);
  int *Message_from_dcminus1nodes = new int[MAXdc];
  int *messages = new int[MAXdc > MAXdv ? MAXdc : MAXdv];
  int *messageFromDvminus1Nodes = new int[MAXdv];
  const char *Files[11]={"data0.tmp","data1.tmp","data2.tmp","data3.tmp",
                  "data4.tmp","data5.tmp","data6.tmp","data7.tmp","data8.tmp","data9.tmp","data10.tmp"};
  FILE *fp = NULL;
  int Nsim=100;
  float pv[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  for(int pe=0;pe<11;pe++){
  float errors[length]={0};
  cout<<pe;
      for(int  num=1;num<=Nsim;num++){
      int *Recieved = BEC_channel(Transmitted, Number_of_VNs, pv[pe]);
      for (int i = 0; i < Number_of_VNs; i++)
      {
        vn = Variable_Nodes + i;
        vn->value = Recieved[i];
        for (int j = 0; j < dv[i]; j++)
          *vn->connected_nodes[j] = vn->value;
      }
      int iteration_index = 0; // Intialisng iteration index to 0
      int max_iterations = 30; // maximum allowed iterations
      while (iteration_index <= max_iterations)
      {
        errors[iteration_index]+=no_of_erasures(Variable_Nodes,Number_of_VNs);
        iteration_index++;
         if(are_Decoded(Variable_Nodes,Number_of_VNs))
         break;
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
      }
      for(int i=0;i<length;i++)
       errors[i]/=Nsim*Number_of_VNs;
      fp=fopen(Files[pe],"w");
      for(int i=0;i<length;i++){
       fprintf(fp,"%d %f\n",i,errors[i]);
       }
      fclose(fp);
      cout<<pe;
  }
   
  FILE *gnupipe = NULL;
  char *GnuCommands[] = {"set title \"Convergence plot for BEC Hard Decision\"",
                          "set key outside",
                         "set yrange [0:1.2]",
                         "set xrange [0:25]",
                         "set xlabel \"iteration index\"",
                         "set ylabel \"error probability\"",
                         "plot 'data0.tmp' with lines title \"p=0.0 \" linecolor 0 linewidth 2,'data1.tmp' with lines title \"p=0.1\" linecolor 1 linewidth 2,'data2.tmp' with lines title \"p=0.2\" linecolor 2 linewidth 2,'data3.tmp' with lines title \"p=0.3\" linecolor 3 linewidth 2,'data4.tmp' with lines title \"p=0.4\" linecolor 4 linewidth 2,'data5.tmp' with lines title \"p=0.5\"linecolor 5 linewidth 2,'data6.tmp' with lines title \"p=0.6\" linecolor 6 linewidth 2,'data7.tmp' with lines title \"p=0.7\" linecolor 7 linewidth 2,'data8.tmp' with lines title \"p=0.8\" linecolor 8 linewidth 2,'data9.tmp' with lines title \"p=0.9\" linecolor 9 linewidth 2,'data10.tmp' with lines title \"p=1.0 \" linecolor 10 linewidth 2",
                         "set terminal wxt size 700,600",
                         "set grid",
                         "replot"};
  gnupipe = _popen("gnuplot -persitent,", "w");
  for (int i = 0; i < sizeof(GnuCommands)/sizeof(GnuCommands[0]); i++)
  {
    fprintf(gnupipe, "%s\n", GnuCommands[i]);
  }

  cout << "Done\n,By this time, the graph should be visible to you";
  return 0;
}
