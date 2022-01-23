#include <bits/stdc++.h>
#include <fstream>
using namespace std;
struct Node
{
    int value;
    double **connected_nodes;
    double *Message_Recieved;
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
bool are_matched(Node *vns, int *Transmitted, int n)
{
    for (int i = 0; i < n; i++)
    {
        if (vns[i].value != Transmitted[i])
            return false;
    }
    return true;
}
double Calculate_CntoVn(double *arr, int n)
{

    double r0giveni = 1;
    for (int i = 0; i < n; i++)
    {
        r0giveni *= (1 - 2 * arr[i]);
    }
    r0giveni = 0.5 + 0.5 * r0giveni;
    return 1 - r0giveni;
}
double Calculate_vntocn(double *arr, int n)
{
    double r = 1;
    for (int i = 0; i < n; i++)
    {
        r *= arr[i];
    }
    return r;
}
int no_of_erasures(Node *Vns,int n)
{
  int ner=0;
  for(int i=0;i<n;i++)
  {
    if(Vns[i].value!=0)
      ner++;
  }
  return ner;
}
int main(int argc, char const *argv[])
{
    ifstream File;
    File.open("HMATRIX_9x12.txt");
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
    for (int i = 0; i < Number_of_CNs; i++)
    {
        Check_Nodes[i].connected_nodes = new double *[dc[i]];
        Check_Nodes[i].Message_Recieved = new double[dc[i]];
        for (int j = 0; j < dc[i]; j++)
            Check_Nodes[i].Message_Recieved[j] = 4;
    }
    for (int i = 0; i < Number_of_VNs; i++)
    {
        Variable_Nodes[i].connected_nodes = new double *[dv[i]];
        Variable_Nodes[i].Message_Recieved = new double[dv[i]];
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
            Node *vn = Variable_Nodes + k;
            int l;
            for (l = 0; l < dv[k]; l++)
            {
                if (vn->Message_Recieved[l] == 4)
                    break;
            }
            Check_Nodes[i].connected_nodes[j] = &vn->Message_Recieved[l];
            *Check_Nodes[i].connected_nodes[j] = 3;
            j++;
            vn = NULL;
            delete vn;
        }
    }
    for (int i = 0; i < Number_of_VNs; i++)
    {
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
    int MAXdc = *max_element(dc, dc + Number_of_CNs);
    int MAXdv = *max_element(dv, dv + Number_of_VNs);
    double *Message_from_dcminus1nodes = new double[MAXdc];
    double *messages = new double[MAXdc > MAXdv ? MAXdc : MAXdv];
    double *messageFromDvminus1Nodes = new double[MAXdv];
    double *pIntial = new double[Number_of_VNs];
    Node *vn = NULL, *cn = NULL;
    int Transmitted[Number_of_VNs] = {0};
    int length = 21;
    
const char *Files[11]={"data0.tmp","data1.tmp","data2.tmp","data3.tmp",
                  "data4.tmp","data5.tmp","data6.tmp","data7.tmp","data8.tmp","data9.tmp","data10.tmp"};
  FILE *fp = NULL;
  int Nsim=100;
  float pv[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  for(int pe=0;pe<11;pe++){
  float errors[length]={0};
  cout<<pe;
    double *q0 = new double[MAXdv];
    double *q1 = new double[MAXdv];
        int Nsim = 1000; // Number of simulations
        for (int num = 1; num <= Nsim; num++)
        {
            int *Recieved = BEC_channel(Transmitted, Number_of_VNs, pv[pe]);
            for (int i = 0; i < Number_of_VNs; i++)
            {
                if (Recieved[i] == 1)
                    pIntial[i] = 1;
                else if (Recieved[i]== -1)
                    pIntial[i] = 0.5;
                else
                    pIntial[i] = 0;
                vn = Variable_Nodes + i;
                vn->value = Recieved[i];
                for (int j = 0; j < dv[i]; j++)
                    *vn->connected_nodes[j] = pIntial[i];
            }
            int iteration_index = 0;
            int Max_iterations = 20;
            while (iteration_index <= Max_iterations)
            {        
                errors[iteration_index]+=no_of_erasures(Variable_Nodes,Number_of_VNs);
                iteration_index++;
                for (int i = 0; i < Number_of_CNs; i++)
                {
                    cn = Check_Nodes + i;
                    for (int j = 0; j < dc[i]; j++)
                        messages[j] = cn->Message_Recieved[j];
                    for (int j = 0; j < dc[i]; j++)
                    {
                        int counter = 0;
                        for (int k = 0; k < dc[i]; k++)
                        {
                            if (k != j)
                            {
                                Message_from_dcminus1nodes[counter++] = messages[k];
                            }
                        }

                        *cn->connected_nodes[j] = Calculate_CntoVn(Message_from_dcminus1nodes, dc[i] - 1);
                    }
                    cn->value = Calculate_CntoVn(messages, dc[i]) >= 0.5 ? 1 : 0;
                }
                bool changed = false;
                for (int i = 0; i < Number_of_VNs; i++)
                {
                    vn = Variable_Nodes + i;
                    for (int j = 0; j < dv[i]; j++)
                        messages[j] = vn->Message_Recieved[j];
                    for (int j = 0; j < dv[i]; j++)
                    {
                        int counter = 0;
                        for (int k = 0; k < dv[i]; k++)
                        {
                            if (k != j)
                            {
                                q0[counter] = 1 - messages[k];
                                q1[counter] = messages[k];
                                counter++;
                            }
                        }
                        double Q0 = (1 - pIntial[i]) * Calculate_vntocn(q0, dv[i] - 1);
                        double Q1 = (pIntial[i]) * Calculate_vntocn(q1, dv[i] - 1);
                        double K = 1 / (Q0 + Q1);
                        *vn->connected_nodes[j] = K * Q1;
                    }
                    for (int j = 0; j < dv[i]; j++)
                    {
                        q0[j] = 1 - messages[j];
                        q1[j] = messages[j];
                    }
                    double Q0 = (1 - pIntial[i]) * Calculate_vntocn(q0, dv[i]);
                    double Q1 = (pIntial[i]) * Calculate_vntocn(q1, dv[i]);
                    int c = vn->value;
                    if(Q1>=Q0) vn->value=1;
                    else if(Q1<Q0) vn->value=0;
                    if (c != vn->value)
                        changed = true;
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
  char *GnuCommands[] = {"set title \"Convergence plot for BEC Soft Decision\"",
                          "set key outside",
                         "set yrange [0:1.2]",
                         "set xrange [0:20]",
                         "set xlabel \"iteration index\"",
                         "set ylabel \"error probability\"",
                         "plot 'data0.tmp' with lines title \"p=0.0 \" linecolor 0 linewidth 2,'data1.tmp' with lines title \"p=0.1\" linecolor 1 linewidth 2,'data2.tmp' with lines title \"p=0.2\" linecolor 2 linewidth 2,'data3.tmp' with lines title \"p=0.3\" linecolor 3 linewidth 2,'data4.tmp' with lines title \"p=0.4\" linecolor 4 linewidth 2,'data5.tmp' with lines title \"p=0.5\"linecolor 5 linewidth 2,'data6.tmp' with lines title \"p=0.6\" linecolor 6 linewidth 2,'data7.tmp' with lines title \"p=0.7\" linecolor 7 linewidth 2,'data8.tmp' with lines title \"p=0.8\" linecolor 8 linewidth 2,'data9.tmp' with lines title \"p=0.9\" linecolor 9 linewidth 2,'data10.tmp' with lines title \"p=1.0 \" linecolor 10 linewidth 2",
                         "set terminal wxt size 700,600",
                         "set grid",
                         "replot"};
  gnupipe = _popen("gnuplot -persitent,", "w");
  for (int i = 0; i < 10; i++)
  {
    fprintf(gnupipe, "%s\n", GnuCommands[i]);
  }
  cout << "Done\n,By this time, the graph should be visible to you";
    return 0;
}
