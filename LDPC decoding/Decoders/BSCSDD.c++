#include <bits/stdc++.h>
using namespace std;
struct Node
{
    int value;
    double **connected_nodes;
    double *Message_Recieved;
};
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
    Node *vn = NULL, *cn = NULL;
    int Transmitted[Number_of_VNs] = {0};
    int length = 11;
    float pVec[length] = {0};
    float pSuccessVec[length] = {0};
    float intial_value = 0, increment = 1 / (float(length - 1));
    double *q0 = new double[MAXdv];
    double *q1 = new double[MAXdv];
    for (int i = 0; i < length; i++)
    {
        pVec[i] = intial_value;
        intial_value += increment;
    }
    double *pIntial=new double[Number_of_VNs];
    for (int pe = 0; pe < length; pe++)
    {
        int Nsim = 1000; // Number of simulations
        int Nerr = 0;
        int Ncorr = 0;
        for (int num = 1; num <= Nsim; num++)
        {
            int *Recieved = Bsc_channel(Transmitted, Number_of_VNs, pVec[pe]);
            for (int i = 0; i < Number_of_VNs; i++)
            {
                if (Recieved[i] == 1)
                    pIntial[i] = 1 - pVec[pe];
                else
                    pIntial[i] = pVec[pe];
                vn = Variable_Nodes + i;
                vn->value = Recieved[i];
                for (int j = 0; j < dv[i]; j++)
                    *vn->connected_nodes[j] = pIntial[i];
            }
            int iteration_index = 0;
            int Max_iterations = 50;
            while (iteration_index <= Max_iterations)
            {
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
                    vn->value = Q1 / Q0 >= 1 ? 1 : 0;
                    if (c != vn->value)
                        changed = true;
                }
                if (!changed)
                    break;
            }

            // After decoding, if it is matchin with transmitted then increase Ncorr
            if (are_matched(Variable_Nodes, Transmitted, Number_of_VNs))
                Ncorr++;
            else
                Nerr++;
        }
        cout << pe;
        pSuccessVec[pe] = (float(Ncorr) / Nsim);
        // storing the value in array for succesfull probability
    }

    // Plotting graph
    FILE *fp = NULL;
    FILE *gnupipe = NULL;
    char *GnuCommands[] = {"set title \"Successful Decoding probability for BSC soft decision\"",
                           "set key noautotitle",
                           "set xlabel \"BSC error probabilty (p)\"",
                           "set ylabel \"successful decoding probability\"",
                           "plot 'data.tmp' with lines linecolor 3 linewidth 2",
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
