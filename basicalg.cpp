#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
using namespace std;
int iAntCount;
int iItCount;
int iLociCount;
int iLociModel;
int iEpiModel;
double alpha;
double rou;
double phe;
double pvalue;
int iTopModel;
int** loci_TopModel;
int iTopLoci;
int* loci_TopLoci;
double* eva_TopModel;
double* phe_TopLoci;
vector<vector<double> > sigepi;
void quickSort(double arr[], int left, int right) {
int i = left;
int j = right;
double tmp;
double pivot = arr[(left + right) / 2];
while (i <= j)
{
while (arr[i] > pivot)
i++;
while (arr[j] < pivot)
j--;
if (i <= j)
{
tmp = arr[i];
arr[i] = arr[j];
arr[j] = tmp;
i++;      
j--;     
}         
}        
if (left < j) 
quickSort(arr, left, j);
if (i < right)
quickSort(arr, i, right);
}
class SNP
{
public:
	int iSample;
	int iLoci;
        vector<int> status;
	vector< vector <int> > data;
	vector<string> SNPnames;
        map<string,int> id;
	int classvalues[2];
        void input_data(char* path);
        void destroy();
        void setpheromone(double level);
        double* pheromone;
        double* cdf;
        void print_inform();
};
void SNP::destroy()
{
delete []pheromone;
delete []cdf;
}
void SNP::setpheromone(double level)
{
int i;
cdf[0]=0;
for(i=0;i<iLoci-1;i++)
{
pheromone[i]=level;
cdf[i+1]=double(i+1)/double(iLoci-1);
}
}
void SNP::print_inform()
{
;
}
void SNP::input_data(char* path)
{
        cout<<"Reading SNP data..."<<endl;
	int i,j,temp;
        string line;
	string word;
        classvalues[0]=0;
        classvalues[1]=0;
	ifstream in1(path);
        i=0;
	while(!in1.eof())
	{
		getline(in1,line);
                if(line=="")
                break;
		istringstream values(line);
                if(i==0)
                {
                   getline(values,word,'\t');
                   while(!values.eof())
                     {
                     getline(values,word,'\t');
                     istringstream int_iss(word);
                     int_iss>>temp;
                     classvalues[temp]++;
                     status.push_back(temp);
                     }
                }
                else
                {
                   getline(values,word,'\t');
                   SNPnames.push_back(word);
                   id[word]=i-1;
                   vector<int>one_locus;
		   while(!values.eof())
		   {
                     getline(values,word,'\t');
                     istringstream int_iss(word);
                     int_iss>>temp;
                     one_locus.push_back(temp);
                   }
                   data.push_back(one_locus);
		}
                i++;
	}
	in1.close();
        iSample=status.size();
        iLoci=SNPnames.size();
        pheromone=new double[iLoci-1];
        cdf=new double[iLoci+1];
        cout<<"-----Reading data completed!------"<<endl;
        cout<<"Number of loci: "<<iLoci<<endl;
        cout<<"Number of samples: "<<iSample<<endl;
}
SNP SNPdata;
class pathway
{
public:
        int num;
        double* pheromone;
        vector<vector <int> >id;
        vector<string> name;
        void input_data(char* path);
        void print_inform();
        void sort_pw(double cutoff);
};
void pathway::print_inform()
{
;
}
void pathway::input_data(char* path)
{
    cout<<"Reading pathway data..."<<endl;
    ifstream in(path);
    string line;
    string word;
    string pwname;
    while(!in.eof())
    {
        getline(in,line);
        if(line=="")
        break;
        istringstream values(line);
        getline(values,pwname,'\t');
        vector<int> one_pw;
        while(!values.eof())
        {
           getline(values,word,'\t');
           map<string,int>::iterator it;
           it=SNPdata.id.find(word);
           if(it!=SNPdata.id.end())
           one_pw.push_back(it->second);
        }
        if(one_pw.size()>1)
        {
        name.push_back(pwname);
        id.push_back(one_pw);
        }
    }
    num=name.size();
}
void pathway::sort_pw(double cutoff)
{
    int i,j,k;
    pheromone=new double[num];
    for(i=0;i<num;i++)
    {
    double* temp;
    j=id[i].size();
    temp=new double[j];
    for(k=0;k<j;k++)
    temp[k]=SNPdata.pheromone[id[i][k]];
    quickSort(temp,0,j-1);
    int endp=(int)((double)j*cutoff);
    if(endp<1)
       endp=1;
    pheromone[i]=0;
    for(k=0;k<endp;k++)
    {
    pheromone[i]+=temp[k];
    }
    pheromone[i]=pheromone[i]/(double)endp;
    delete[] temp;
    }
    for(i=num-1;i>0;i--)
    {
       for(j=0;j<i;j++)
       {
       if(pheromone[j]<pheromone[j+1])
          {
          double tp=pheromone[j];
          string tn=name[j];
          pheromone[j]=pheromone[j+1];
          name[j]=name[j+1];
          pheromone[j+1]=tp;
          name[j+1]=tn;
          }
       }
    }
}
pathway pdata;
