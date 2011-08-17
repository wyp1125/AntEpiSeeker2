#include "basicalg.cpp"
double chi_square(int* selectedSNPSet,int k)
{
int comb = (int)pow(3.0, k);
double** observedValues;
double* colSumTable;
double** expectedValues;
int i,j,index;
observedValues=new double*[2];
expectedValues=new double*[2];
for(i=0;i<2;i++)
{
observedValues[i]=new double[comb];
expectedValues[i]=new double[comb];
}
colSumTable=new double[comb];
for(i=0;i<comb;i++)
{
observedValues[0][i] = 0;
observedValues[1][i] = 0;
colSumTable[i] = 0;
}
/*constructing observed freq table*/
bool cont;
for(i=0;i<SNPdata.iSample;i++)
{
index = 0;
cont = 1;
for(j=0;j<k;j++)
{
if(SNPdata.data[selectedSNPSet[j]][i]==3)
{
cont = 0;
break;
}
else
{
index+=SNPdata.data[selectedSNPSet[j]][i]*(int)pow(3.0,(k-1-j));
}
}
if(cont){
if(!SNPdata.status[i]==0)
{						
observedValues[0][index]++;
}
else{
observedValues[1][index]++;
}
colSumTable[index]++;
}
}
/*computing expected freq values and compute chi-square value*/
double x2 = 0;
for(i=0;i<comb;i++){
expectedValues[0][i] = colSumTable[i]*SNPdata.classvalues[0]/(double)SNPdata.iSample;
expectedValues[1][i] = colSumTable[i]*SNPdata.classvalues[1]/(double)SNPdata.iSample;
if(expectedValues[0][i]!=0){
x2 = x2 + (expectedValues[0][i]-observedValues[0][i])*(expectedValues[0][i]-observedValues[0][i])/expectedValues[0][i];
}
if(expectedValues[1][i]!=0){
x2 = x2 + (expectedValues[1][i]-observedValues[1][i])*(expectedValues[1][i]-observedValues[1][i])/expectedValues[1][i];
}
}
for(i=0;i<2;i++)
{
delete [] expectedValues[i];
delete [] observedValues[i];
}
delete [] expectedValues;
delete [] observedValues;
delete [] colSumTable;
return x2; 
}
