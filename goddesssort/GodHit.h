using namespace std;
#define N 50
class GSHit
{
  public:
  int nHits; 
  int ID[N];
  double E[N];
  double T[N];
  double X[N];
  double Y[N];
  double Z[N];
  GSHit(){
    int i;
    for (i=0;i<N;i++) {
      nHits=0;ID[i]=0;E[i]=0;T[i]=0;X[i]=0;Y[i]=0;Z[i]=0;
    }
  };
  void Reset()
  {
    int i;
    for (i=0;i<N;i++) {
      nHits=0;ID[i]=0;E[i]=0;T[i]=0;X[i]=0;Y[i]=0;Z[i]=0;
    }
  };
};

class BGOHit
{
  public:
  int nHits; 
  int ID[N];
  double E[N];
  double T[N];
  double X[N];
  double Y[N];
  double Z[N];
  BGOHit(){
    int i;
    for (i=0;i<N;i++) {
      nHits=0;ID[i]=0;E[i]=0;T[i]=0;X[i]=0;Y[i]=0;Z[i]=0;
    }
  };
  void Reset()
  {
    int i;
    for (i=0;i<N;i++) {
      nHits=0;ID[i]=0;E[i]=0;T[i]=0;X[i]=0;Y[i]=0;Z[i]=0;
    }
  };
};

class GodEvent
{
  public:
  GodEvent(){};
  int nHits;
  GSHit theGS;
  BGOHit theBGO;
  void Reset()
  {
    theGS.Reset();
    theBGO.Reset();
  };
};


