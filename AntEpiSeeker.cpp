#include <gsl/gsl_cdf.h>
#include <ctime>
#include <cstdlib>
#include "model.cpp"
int rnd(int uper)
{
          return (rand()%uper);
}
double rnd(int low, double uper)
{
double p=(rand()/(double)RAND_MAX)*((uper)-(low))+(low);
return (p);
}
int cdf2locus(double x, int start, int end)
{
if(start+1==end)
{
return start;
}
else
{
int temp=(start+end)/2;
if(SNPdata.cdf[temp]<=x)
{
return cdf2locus(x,temp,end);
}
else
{
return cdf2locus(x,start,temp);
}
}
}
class ant
{
private:
 int m_iLociCount;
public:
 int ant_number;
 int* tabu;
 void initiate();
 void start();
 void move();
 int ChooseNextLocus();
 void addlocus(int locus);
 void destroy();
};
void ant::initiate()
{
tabu=new int[iLociModel+1];
}
void ant::start()
{
 m_iLociCount=0;
}
int ant::ChooseNextLocus()
{
int i,j,k;
j=-1;
double mRate;
while(j==-1)
{
j=0;
srand((unsigned)time(NULL)+rand());
mRate=rnd(0,1);
k=cdf2locus(mRate,0,iLociCount-1);
for(i=0;i<m_iLociCount;i++)
{
if(k==tabu[i])
{
j=-1;
break;
}
}
}
return k;
}
void ant::addlocus(int locus)
{
 tabu[m_iLociCount]=locus;
 m_iLociCount++;
}
void ant::move()
{
int j;
j=ChooseNextLocus();
addlocus(j);
}
void ant::destroy()
{
delete[] tabu;
}
class project
{
public:
project();
~project();
ant* ants;
void GetAnt();
void UpdatePheromone();
void StartSearch();
};
project::project()
{
ants=new ant[iAntCount];
int i;
for(i=0;i<iAntCount;i++)
{ants[i].initiate();}
}
project::~project()
{
delete [] ants;
}
void project::UpdatePheromone()
{
int i,j,k,tag;
for(i=0;i<iLociCount;i++)
SNPdata.pheromone[i]=SNPdata.pheromone[i]*(1-rou);
double eva;
int* locidata;
locidata=new int[iLociModel];
for(i=0;i<iAntCount;i++)
{
for(j=1;j<=iLociModel;j++)
{
locidata[j-1]=ants[i].tabu[j];
}
eva=chi_square(locidata,iLociModel)/100;
for(j=0;j<iLociModel;j++)
{
SNPdata.pheromone[locidata[j]]+=eva;
}
if(eva>eva_TopModel[0])
{
tag=1;
for(j=0;j<iTopModel;j++)
{
if(fabs(eva-eva_TopModel[j])<0.000001)
{
tag=0;
break;
}
}
if(tag)
{
eva_TopModel[0]=eva_TopModel[1];
for(j=0;j<iLociModel;j++)
loci_TopModel[0][j]=loci_TopModel[1][j];
k=1;
while(k<iTopModel && eva>eva_TopModel[k])
{
eva_TopModel[k-1]=eva_TopModel[k];
for(j=0;j<iLociModel;j++)
loci_TopModel[k-1][j]=loci_TopModel[k][j];
k++;
}
eva_TopModel[k-1]=eva;
for(j=0;j<iLociModel;j++)
loci_TopModel[k-1][j]=locidata[j];
}
}
}
//#################update cdf##################
double * prob;
double temp=0;
prob=new double[iLociCount];
for(i=0;i<iLociCount;i++)
{
temp+=pow(SNPdata.pheromone[i],alpha);
prob[i]=temp;
}
for(i=0;i<iLociCount-1;i++)
{
SNPdata.cdf[i+1]=prob[i]/temp;
}
delete [] prob;
delete [] locidata;
}
void project::GetAnt()
{
 int i=0;
 int locus;
 srand((unsigned)time(NULL)+rand());
 for (i=0;i<iAntCount;i++)
 {
  ants[i].start();
  locus=rnd(iLociCount);
  ants[i].ant_number=i;
  ants[i].addlocus(locus);
 }
 }
void project::StartSearch()
{
int max,i,j;
max=0;
while(max<iItCount)
{
cout<<"Iteration: "<<max<<endl;
GetAnt();
for(i=0;i<iAntCount;i++)
{
for(j=1;j<=iLociModel;j++)
{
ants[i].move();
}
}
UpdatePheromone();
max++;
}
}
void get_toploci()
{
int i,j,temp1;
double temp2;
int* tag;
tag=new int[iLociCount];
double* SortPheromone;
SortPheromone=new double[iLociCount];
for(i=0;i<iLociCount;i++)
{
tag[i]=i;
SortPheromone[i]=SNPdata.pheromone[i];
}
for(i=0;i<iTopLoci;i++)
{
for(j=i+1;j<iLociCount;j++)
{
if(SortPheromone[j]>SortPheromone[i])
{
temp1=tag[i];
tag[i]=tag[j];
tag[j]=temp1;
temp2=SortPheromone[i];
SortPheromone[i]=SortPheromone[j];
SortPheromone[j]=temp2;
}
}
loci_TopLoci[i]=tag[i];
phe_TopLoci[i]=SortPheromone[i];
}
delete []tag;
delete []SortPheromone;
}
char pinputfile[200];
char inputfile[200];
char outputfile[200];
char pwyoutfile[200];
double rnd_values[10000];
int ihapsize[2];
int counts;
double pcut;
vector<vector< double> > interactions;
vector<vector< double> > mini_interactions;
void postprocessing();
void mini_fp();
void write_result(char* path,int mod);
void loadparameters(char* path);
int main()
{
        int i,a;
        iAntCount=1000;
        counts=300;
        ihapsize[0]=6;
        ihapsize[1]=3;
        alpha=1;
        rou=0.05;
        phe=100;
        iTopModel=1000;
        iTopLoci=200;
        iEpiModel=2;
        pvalue=0.01;
        char* para="parameters.txt";
        loadparameters(para);
        SNPdata.input_data(inputfile);
        SNPdata.setpheromone(phe);
        pdata.input_data(pinputfile);
        iLociCount=SNPdata.iLoci;
        char* logpath="AntEpiSeeker.log";
        char* minipath="results_maximized.txt";
        write_result(logpath,0);
        for(a=ihapsize[0];a>=ihapsize[1];a--)
        {
        iItCount=counts;
        cout<<ihapsize[0]-a+1<<" round search"<<endl;
        iLociModel=a;
        loci_TopModel=new int* [iTopModel];
        for(i=0;i<iTopModel;i++)
                loci_TopModel[i]=new int [iLociModel];
        loci_TopLoci=new int [iTopLoci];
        phe_TopLoci=new double [iTopLoci];
        eva_TopModel=new double [iTopModel];
        for(i=0;i<iTopModel;i++)
            eva_TopModel[i]=0;
        cout<<"-----Initializing parameters-----"<<endl;
        cout<<"Number of ants: "<<iAntCount<<endl;
        cout<<"Number of iterations: "<<iItCount<<endl;
        cout<<"Start pheromone level: "<<phe<<endl;
        cout<<"Rou: "<<rou<<endl;
        cout<<"Alpha: "<<alpha<<endl;
        cout<<"Number of top ranking SNP sets: "<<iTopModel<<endl;
        cout<<"Number of top ranking loci: "<<iTopLoci<<endl;
        cout<<"Size of SNP sets: "<<a<<endl;
        cout<<"Number of SNPs in an epistatic interaction: "<<iEpiModel<<endl;
        project episeeker;
        episeeker.StartSearch();
        get_toploci();
        write_result(logpath,3);
        postprocessing();
        for(i=0;i<iTopModel;i++)
        {
         delete [] loci_TopModel[i];
        }
         delete [] loci_TopModel;
         delete [] loci_TopLoci;
         delete [] phe_TopLoci;
         delete [] eva_TopModel;
         
        }
        mini_fp();
        write_result(minipath,4);
        write_result(outputfile,1);
        cout<<"Sorting pathways..."<<endl;
        pdata.sort_pw(pcut);
        cout<<"Done!"<<endl;
        write_result(pwyoutfile,2);
        SNPdata.destroy();
        return 1;
}
void mini_fp()
{
int i,j,k,l,find,hit,first_pos;
if(interactions.size()>0)
{
mini_interactions.push_back(interactions[0]);
for(i=1;i<interactions.size();i++)
{
find=0;
for(j=0;j<mini_interactions.size();j++)
{
hit=0;
for(k=0;k<iEpiModel;k++)
{
for(l=0;l<iEpiModel;l++)
{
if(interactions[i][k]==mini_interactions[j][l])
{
hit=1;
find++;
break;
}
}
if(hit==1)
break;
}
if(find==1 && hit==1)
{
if(interactions[i][iEpiModel]<=mini_interactions[j][iEpiModel])
{
break;
}
else
{
for(k=0;k<iEpiModel+2;k++)
{
mini_interactions[j][k]=interactions[i][k];
}
first_pos=j;
}
}
if(find>1 && hit==1)
{
if(mini_interactions[j][iEpiModel]>mini_interactions[first_pos][iEpiModel])
{
mini_interactions.erase(mini_interactions.begin()+first_pos);
break;
}
else
{
mini_interactions.erase(mini_interactions.begin()+j);
break;
}
}
}
if(find==0)
{
mini_interactions.push_back(interactions[i]);
}
}
///////////////////////////////////////////////////////////////
}
}
void loadparameters(char* path)
{
        FILE *f = fopen(path, "r");
        if(f == NULL)
        {       printf("Cannot open the parameter file \"%s\"\n", path);
                return;
        }
        int i;
        char tmp[1000];
        while(fgets(tmp, 1000, f) != NULL)
        {
                char str[1000];
                sprintf(str, "%s", tmp);
                str[13]=0;
                if(strcmp(str, "iItCountHsize") == 0)
                {
                        for(i = 13; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        counts=atoi(&tmp[i]);
                }
                str[12]=0;
                if(strcmp(str, "largesetsize") == 0)
                {
                        for(i = 12; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        ihapsize[0]=atoi(&tmp[i]);
                }
                if(strcmp(str, "smallsetsize") == 0)
                {
                        for(i = 12; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        ihapsize[1]=atoi(&tmp[i]);
                }
                str[9] = 0;
                if(strcmp(str, "iTopModel") == 0)
                {       for(i = 9; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        iTopModel=atoi(&tmp[i]);
                }
                else if(strcmp(str, "iEpiModel") == 0)
                {       for(i = 9; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        iEpiModel=atoi(&tmp[i]);
                }
                else if(strcmp(str, "iAntCount") == 0)
                {       for(i = 9; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        iAntCount=atoi(&tmp[i]);
                }
                str[8]=0;
                if(strcmp(str, "iTopLoci") == 0)
                {       for(i = 8; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        iTopLoci=atoi(&tmp[i]);
                }
                str[7] = 0;
                if(strcmp(str, "INPFILE") == 0)
                {
                for(i = 7; i < (int)strlen(tmp); i++)
                    if(tmp[i] != ' ' && tmp[i] != '\t')
                            break;
                    int j;
                    for(j = i + 1; j < (int)strlen(tmp); j++)
                         if(tmp[j] == ' ' || tmp[j] == '\t' || tmp[j] == '"' || tmp[j]=='\n')
                         break;
                    tmp[j] = 0;
                    if(tmp[i] == '"') sprintf(inputfile, "%s", &tmp[i + 1]);
                        else sprintf(inputfile, "%s", &tmp[i]);
                }
                if(strcmp(str, "OUTFILE") == 0)
                {
                for(i = 7; i < (int)strlen(tmp); i++)
                    if(tmp[i] != ' ' && tmp[i] != '\t')
                         break;
                    int j;
                    for(j = i + 1; j < (int)strlen(tmp); j++)
                        if(tmp[j] == ' ' || tmp[j] == '\t' || tmp[j] == '"' || tmp[j]=='\n')
                        break;
                    tmp[j] = 0;
                    if(tmp[i] == '"') sprintf(outputfile, "%s", &tmp[i + 1]);
                        else sprintf(outputfile, "%s", &tmp[i]);
                }
                if(strcmp(str, "PWSNPFL") == 0)
                {
                for(i = 7; i < (int)strlen(tmp); i++)
                    if(tmp[i] != ' ' && tmp[i] != '\t')
                        break;
                    int j;
                    for(j = i + 1; j < (int)strlen(tmp); j++)
                        if(tmp[j] == ' ' || tmp[j] == '\t' || tmp[j] == '"' || tmp[j]=='\n')
                        break;
                    tmp[j] = 0;
                    if(tmp[i] == '"') sprintf(pinputfile, "%s", &tmp[i + 1]);
                        else sprintf(pinputfile, "%s", &tmp[i]);
                }
                if(strcmp(str, "PWYFILE") == 0)
                {
                for(i = 7; i < (int)strlen(tmp); i++)                                              
                    if(tmp[i] != ' ' && tmp[i] != '\t')
                        break;
                    int j;
                    for(j = i + 1; j < (int)strlen(tmp); j++)
                        if(tmp[j] == ' ' || tmp[j] == '\t' || tmp[j] == '"' || tmp[j]=='\n')                                        
                        break;
                    tmp[j] = 0;
                    if(tmp[i] == '"') sprintf(pwyoutfile, "%s", &tmp[i + 1]);
                        else sprintf(pwyoutfile, "%s", &tmp[i]);
                }
                str[6]=0;
                if(strcmp(str, "pvalue") == 0)
                {       for(i = 6; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        pvalue = atof(&tmp[i]);
                }
                if(strcmp(str, "pwprop") == 0)
                {       for(i = 6; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        pcut = atof(&tmp[i]);
                }
                str[5] = 0;
                if(strcmp(str, "alpha") == 0)
                {
                        for(i = 5; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        alpha=atof(&tmp[i]);
                }
                str[3]=0;
                if(strcmp(str, "rou") == 0)
                {       for(i = 3; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        rou = atof(&tmp[i]);
                }
                else if(strcmp(str, "tau") == 0)
                {       for(i = 3; i < (int)strlen(tmp); i++)
                                if(tmp[i] >= 48 && tmp[i] < 58)
                                        break;
                        phe = atof(&tmp[i]);
                }
        }
        fclose(f);
 
}
int **throughout(int m,int n,int** result)
{
int num,mi,ni,i,j,temp;
num=0;
mi=m-1;
ni=n-1;
int* comb=new int[n];
for(i=0;i<=ni;i++)
{
comb[i]=i;
}
temp=0;
while(1)
{
temp++;
while(comb[ni]==mi)
{
mi--;
ni--;
if(ni==-1)
{
break;
}
}
if(ni==-1)
{
break;
}
if(ni<n-1)
{
comb[ni]++;
for(j=ni+1;j<n;j++)
{
comb[j]=comb[j-1]+1;
}
mi=m-1;
ni=n-1;
}
while(comb[ni]<m)
{
for(i=0;i<n;i++)
{
result[num][i]=comb[i];
}
num++;
comb[ni]++;
}
comb[ni]--;
}
delete comb;
return result;
}
long int comb_num(int m,int n)
{
long int i,p,q;
p=1;
q=1;
for(i=1;i<=n;i++)
{
p=p*i;
q=q*(m-i+1);
}
return q/p;
}  
void postprocessing()
{
long int i,num;
int j,k,l,s,isNew;
double eva,p_value;
int* loci;
loci=new int[iEpiModel];
vector<double> record;
record.resize(iEpiModel+2);
cout<<"Post-processing"<<endl;
cout<<"Dealing with top ranking SNP sets"<<endl;
num=comb_num(iLociModel,iEpiModel);
cout<<"Number of interatcions evaluated: "<<num*iTopModel<<endl;
int** combination;
combination=new int*[num];
for(i=0;i<num;i++)
combination[i]=new int[iEpiModel];
combination=throughout(iLociModel,iEpiModel,combination);
for(i=iTopModel-1;i>=0;i--)
{
for(j=0;j<num;j++)
{
for(k=0;k<iEpiModel;k++)
{
loci[k]=loci_TopModel[i][combination[j][k]];
}
eva=chi_square(loci,iEpiModel);
p_value=1-gsl_cdf_chisq_P(eva,pow(3,iEpiModel)-1);
if(p_value<pvalue)
{
isNew=1;
for(l=0;l<interactions.size();l++)
{
if(fabs(eva-interactions[l][iEpiModel])<0.00000001)
{
isNew=0;
break;
}
}
if(isNew)
{
for(s=0;s<iEpiModel;s++)
{
record[s]=(double)loci[s];
}
record[iEpiModel]=eva;
record[iEpiModel+1]=p_value;
interactions.push_back(record);
}
}
}
}
for(i=0;i<num;i++)
delete []combination[i];
delete []combination;
cout<<"Dealing with top ranking loci"<<endl;
num=comb_num(iTopLoci,iEpiModel);
int** combination1;
combination1=new int*[num];
for(i=0;i<num;i++)
combination1[i]=new int[iEpiModel];
combination1=throughout(iTopLoci,iEpiModel,combination1);
cout<<"Number of interactions evaluated: "<<num<<endl;
for(i=0;i<num;i++)
{
for(j=0;j<iEpiModel;j++)
{
loci[j]=loci_TopLoci[combination1[i][j]];
}
eva=chi_square(loci,iEpiModel); 
p_value=1-gsl_cdf_chisq_P(eva,pow(3,iEpiModel)-1); 
if(p_value<pvalue)
{
isNew=1;
for(l=0;l<interactions.size();l++)
{
if(fabs(eva-interactions[l][iEpiModel])<0.00000001)
{
isNew=0;
break;
}
}
if(isNew)
{
for(s=0;s<iEpiModel;s++)
{
record[s]=(double)loci[s];  
}
record[iEpiModel]=eva; 
record[iEpiModel+1]=p_value;  
interactions.push_back(record);
}
}
}
for(i=0;i<num;i++)
delete []combination1[i];
delete []combination1;
}
void write_result(char* path,int mod)
{
int i,j;
std::ofstream result;
if(mod==0)
{
result.open(path,ios::out);
}
if(mod==1)
{
result.open(path,ios::out);
result<<"Epistatic interactions:"<<endl;
result<<"Loci\tChi-square\tP value"<<endl;
for(i=0;i<mini_interactions.size();i++)
{
/////////////assigning p values////////////////
    for(j=0;j<iEpiModel;j++)
    {
        result<<mini_interactions[i][j]<<"("<<SNPdata.SNPnames[int(mini_interactions[i][j])]<<")"<<" ";
    }
    result<<"\t"<<mini_interactions[i][iEpiModel]<<"\t"<<mini_interactions[i][iEpiModel+1]<<endl;         
}
}
if(mod==2)
{
result.open(path,ios::out);
result<<"Rank\tPathway\tPheromone"<<endl;
for(i=0;i<pdata.num;i++)
{
result<<i+1<<"\t"<<pdata.name[i]<<"\t"<<pdata.pheromone[i]<<endl;
}
}
if(mod==3)
{
result.open(path,ios::app);
result<<"#####################Intermediate results############################"<<endl;
result<<"Top ranking SNP sets"<<endl;
for(i=iTopModel-1;i>=0;i--)
{
    for(j=0;j<iLociModel;j++)
    result<<loci_TopModel[i][j]<<" ";
    result<<eva_TopModel[i]<<endl;
}
result<<"Top ranking loci"<<endl;
for(i=0;i<iTopLoci;i++)
{
    result<<SNPdata.SNPnames[loci_TopLoci[i]]<<" "<<phe_TopLoci[i]<<endl;
}
}
if(mod==4)
{
result.open(path,ios::out);
result<<"Epistatic interactions:"<<endl;
result<<"Loci\tChi-square\tP value"<<endl;
for(i=0;i<interactions.size();i++)
{
   for(j=0;j<iEpiModel;j++)
   {
       result<<interactions[i][j]<<"("<<SNPdata.SNPnames[int(interactions[i][j])]<<")"<<" ";
   }
   result<<"\t"<<interactions[i][iEpiModel]<<"\t"<<interactions[i][iEpiModel+1]<<endl;
}
}
result.close();
}
