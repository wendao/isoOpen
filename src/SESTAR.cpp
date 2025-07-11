#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <cmath>
#include <thread>
#include <mutex>
#include <random>
#include <sstream>
#include <iomanip>
#include <stdexcept>

using namespace std;
typedef unsigned char BYTE;

#define PI 3.14159
#define ISOLENGTH 15
#define NEUTRON 1.0033548
#define HPLUS 1.0072765
#define HNEUTRON 1.006277
#define SECUTOFF 10.0
#define RATIOCUTOFF 6.0
#define ENVELOPEPPM 0.000010
#define ELEMENTISOLENGTH 10

const int STATELENGTH = sizeof(int) + sizeof(unsigned short);
const int DATALENGTH = 2 * sizeof(float) + 15 * sizeof(int) + 2 * (15 * sizeof(float));
const int REPORTLENGTH = STATELENGTH + 2 * DATALENGTH + sizeof(bool);
const int RECORD_SIZE = REPORTLENGTH + sizeof(bool);

const int EC = 6; //number of elements using in the distribution calculation
unsigned short* AcceptCharge = new unsigned short[6]{ 6,5,4,3,2,1 };
int AcceptChargeLength = 6;

int ElementNumber = 5;

//C,H,O,N,S
float (*Element)[ELEMENTISOLENGTH] = new float[ElementNumber][ELEMENTISOLENGTH] {
    { 0.989,0.011,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 },
    { 0.99985,0.00015,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 },
    { 0.998,0.0,0.002,0.0,0.0,0.0,0.0,0.0,0.0,0.0 },
    { 0.9963,0.0037,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 },
    { 0.9579,0.0,0.0421,0.0,0.0,0.0,0.0,0.0,0.0,0.0 },
};

double* ExactMass = NULL;
float* InitialCutoff = NULL;

enum ElementName {
	Se=0, Dimethyl=1, Dimethyl2=2, Dimethyl3=3
};

struct EnvelopeParameters {
	int Length;
	int Start;
	int Mass;
	bool operator<(const EnvelopeParameters& a) const {
		return Mass < a.Mass ? true : (Mass > a.Mass ? false : (Length < a.Length ? true : (Length > a.Length ? false : (Start < a.Start))));
	}
};

map<EnvelopeParameters, float> CutoffMap;
vector<set<EnvelopeParameters>> AcceptedParameters;

void AddElement(float distribution[ELEMENTISOLENGTH], double exactmass, float initialcutoff) {
    float(*temp)[ELEMENTISOLENGTH] = new float[ElementNumber + 1][ELEMENTISOLENGTH];
    memcpy(temp, Element, ElementNumber * ELEMENTISOLENGTH * sizeof(float));
    memcpy(temp + ElementNumber, distribution, ELEMENTISOLENGTH * sizeof(float));

    delete[] Element;  // Free old memory
    Element = temp;

    double* ttmp = new double[ElementNumber - 4];
    if (ExactMass) {  // Check if ExactMass exists before copying
        memcpy(ttmp, ExactMass, (ElementNumber-5) * sizeof(double));
        delete[] ExactMass;  // Free old memory
    }
    ttmp[ElementNumber-5] = exactmass;
    ExactMass = ttmp;

    float* tttmp = new float[ElementNumber - 4];
    if (InitialCutoff) {  // Check if InitialCutoff exists before copying
        memcpy(tttmp, InitialCutoff, (ElementNumber-5) * sizeof(float));
        delete[] InitialCutoff;  // Free old memory
    }
    tttmp[ElementNumber-5] = initialcutoff;
    InitialCutoff = tttmp;

    ++ElementNumber;
    AcceptedParameters.push_back(set<EnvelopeParameters>());
}

struct Rscore {
	float R;
	int Position;
	float IsotopeDistribution[ISOLENGTH];
};

class Peptide {
public:
	Peptide(float mass, int specialElementNumber = 0)
	{
		static const float averaginemass = 111.1254;
		static const float composition[EC] = { 4.9348,7.7583,1.4773,1.3577,0.0,0.0 };
		float x = mass / averaginemass;
		for (int i = 0; i < EC; i++) {
			ElementCount[i] = round(composition[i] * x);
		}
		ElementCount[5] = specialElementNumber;
	}

	void CalcDistribution(float IsotopeDistribution[ISOLENGTH], int element) const {
		for (int i = 0; i < ISOLENGTH; i++) {
			IsotopeDistribution[i] = 0;
		}
		IsotopeDistribution[0] = 1;

		for (int i = 0; i < EC - 1; i++) {
			float tmp[ISOLENGTH]{ 0 };
			memcpy(tmp, Element[i], ELEMENTISOLENGTH * sizeof(float));
			RepeatConvolve(tmp, ISOLENGTH, ElementCount[i]);
			Convolve(IsotopeDistribution, ISOLENGTH, tmp, ISOLENGTH);
		}
		float tmp[ISOLENGTH]{ 0 };
		memcpy(tmp, Element[EC - 1 + element], ELEMENTISOLENGTH * sizeof(float));
		RepeatConvolve(tmp, ISOLENGTH, ElementCount[EC - 1]);
		Convolve(IsotopeDistribution, ISOLENGTH, tmp, ISOLENGTH);
	}

private:
	int ElementCount[EC];
	void Convolve(float a[], int lena, float b[], int lenb) const {
		float* r = new float[lena] { 0 };
		for (int i = 0; i < lena; i++) {
			for (int j = 0; j < lenb; j++)
			{
				if (i + j < lena)
				{
					r[i + j] += a[i] * b[j];
				}
				else { break; }
			}
		}
		memcpy(a, r, lena * sizeof(float));
		delete[] r;
		r = NULL;
	}

	void RepeatConvolve(float a[], int length, unsigned int n) const {
		float* tmp = new float[length] { 0 };
		tmp[0] = 1;
		while (n) {
			if (n & 1)
				Convolve(tmp, length, a, length);
			Convolve(a, length, a, length);
			n >>= 1;
		}
		memcpy(a, tmp, length * sizeof(float));
		delete[] tmp;
		tmp = NULL;
	}
};


class Envelope {
public:
	int Tag = 0b10000000000000000000000000000000;

	float Error = -1;
	Envelope* Pair = NULL;

	Envelope(double mz, float intensity , unsigned int position , unsigned short charge) {
		Add(mz, intensity , position);
		Charge = charge;
	}

	void Add(double mz, float intensity , unsigned int position) {
		if (Length >= MaxLength) {
			int last = MaxLength;
			MaxLength += ISOLENGTH;
			double* mtemp = new double[MaxLength] { 0 };
			float* itemp = new float[MaxLength] { 0 };
			int* ptemp = new int[MaxLength] { 0 };
			memcpy(ptemp, Position, last * sizeof(int));
			for (int i = MaxLength - ISOLENGTH; i < MaxLength; i++) {
				ptemp[i] = -1;
			}
			delete[] Position;
			Position = ptemp;
			ptemp = NULL;
			memcpy(mtemp, Mz, last * sizeof(double));
			memcpy(itemp, Intensity, last * sizeof(float));
			
			delete[] Mz;
			delete[] Intensity;
			
			Mz = mtemp;
			Intensity = itemp;
			
			itemp = NULL;
			mtemp = NULL;
		}
		Mz[Length] = mz;
		Intensity[Length] = intensity;
		Position[Length] = position;
		Length += 1;
	}

	BYTE* Report() const {
        	BYTE* r = new BYTE[REPORTLENGTH];
		BYTE* temp = r;
		memcpy(temp, &Tag, sizeof(Tag));
		temp += sizeof(Tag);
		memcpy(temp, &Charge, sizeof(Charge));
		temp += sizeof(Charge);

		temp = DataCopy(temp);

		if (Pair != NULL) {
			*temp = (BYTE)true;
			temp += sizeof(bool);
			temp = Pair->DataCopy(temp);
		}
		else *temp = (BYTE)false;
		temp = NULL;
		
		return r;
	}


	BYTE* Report(int position) const {
        	BYTE* r = new BYTE[REPORTLENGTH];

		BYTE* temp = r;
		memcpy(temp, &Tag, sizeof(Tag));
		temp += sizeof(Tag);
		memcpy(temp, &Charge, sizeof(Charge));
		temp += sizeof(Charge);

		temp = DataCopy(temp);

		if (Pair != NULL) {
			*temp = (BYTE)true;
			temp += sizeof(bool);
			temp = Pair->DataCopy(temp);
		}
		else *temp = (BYTE)false;
		temp = NULL;

		return r;
	}

	inline BYTE* DataCopy(BYTE* ptr) const {
		static const int size = ISOLENGTH * sizeof(float);

		BYTE* temp = ptr;
		float tempf;
		tempf = R[0].R;
		memcpy(temp, &tempf, sizeof(float));
		temp += sizeof(float);
		tempf = R[1].R;
		memcpy(temp, &tempf, sizeof(float));
		temp += sizeof(float);
		memcpy(temp, Position, sizeof(*Position) * ISOLENGTH);
		temp += sizeof(*Position) * ISOLENGTH;
		
		for (int i = 0; i < ISOLENGTH; ++i) {
			tempf = R[0].IsotopeDistribution[i];
			memcpy(temp, &tempf, sizeof(float));
			tempf = R[1].IsotopeDistribution[i];
			memcpy(temp + size, &tempf, sizeof(float));
			temp += sizeof(float);
		}

		return temp + size;
	}

	int GetPosition(int iter) const {
		if (iter >= Length) return -1;
		else return Position[iter];
	}

    // Add proper copy assignment operator
    Envelope& operator=(const Envelope& e) {
        if (this != &e) {
            delete[] Mz;
            delete[] Intensity;
            delete[] Position;
            
            MaxLength = e.MaxLength;
            Length = e.Length;
            Charge = e.Charge;

            Mz = new double[MaxLength];
            Intensity = new float[MaxLength];
            Position = new int[MaxLength];
            
            memcpy(Mz, e.Mz, MaxLength * sizeof(double));
            memcpy(Intensity, e.Intensity, MaxLength * sizeof(float));
            memcpy(Position, e.Position, MaxLength * sizeof(int));
            
            R[0] = e.R[0];
            R[1] = e.R[1];
            Tag = e.Tag;
            Error = e.Error;
            Pair = e.Pair;
        }
        return *this;
    }
    
    // Add move constructor
    Envelope(Envelope&& e) noexcept {
        Mz = e.Mz;
        Intensity = e.Intensity;
        Position = e.Position;
        MaxLength = e.MaxLength;
        Length = e.Length;
        Charge = e.Charge;
        R[0] = e.R[0];
        R[1] = e.R[1];
        Tag = e.Tag;
        Error = e.Error;
        Pair = e.Pair;
        
        e.Mz = nullptr;
        e.Intensity = nullptr;
        e.Position = nullptr;
    }
    
    // Add move assignment operator
    Envelope& operator=(Envelope&& e) noexcept {
        if (this != &e) {
            delete[] Mz;
            delete[] Intensity;
            delete[] Position;
            
            Mz = e.Mz;
            Intensity = e.Intensity;
            Position = e.Position;
            MaxLength = e.MaxLength;
            Length = e.Length;
            Charge = e.Charge;
            R[0] = e.R[0];
            R[1] = e.R[1];
            Tag = e.Tag;
            Error = e.Error;
            Pair = e.Pair;
            
            e.Mz = nullptr;
            e.Intensity = nullptr;
            e.Position = nullptr;
        }
        return *this;
    }

	~Envelope() {
		delete[] Mz;
		Mz = NULL;
		delete[] Intensity;
		Intensity = NULL;
		delete[] Position;
		Position = NULL;
	}

	inline double GetMz(int iter) const {
		return Mz[iter];
	}

	inline float GetIntensity(int iter) const {
		return Intensity[iter];
	}

	inline double GetLastMz() const {
		return Mz[Length - 1];
	}

	inline int GetLength() const {
		return Length;
	}

	inline int GetCharge() const {
		return Charge;
	}

	inline int GetMaxIntensityPosition() const {
		return max_element(Intensity, Intensity + Length) - Intensity;
	}

	inline float GetMaxIntensity() const {
		return *max_element(Intensity, Intensity + Length);
	}

	inline double GetMaxMz() const {
		return Mz[GetMaxIntensityPosition()];
	}

	inline Rscore GetR(int iter) const {
		if (iter > 1) return Rscore{ -1, -1 };
		return R[iter];
	}

	Envelope(const Envelope& e) {
		MaxLength = e.MaxLength;
		Length = e.Length;
		Charge = e.Charge;

		Mz = new double[MaxLength];
		Intensity = new float[MaxLength];
		memcpy(Mz, e.Mz, MaxLength * sizeof(double));
		memcpy(Intensity, e.Intensity, MaxLength * sizeof(float));

		Position = new int[MaxLength];
		memcpy(Position, e.Position, MaxLength * sizeof(int));
	}

	Envelope(const Envelope& e, int start, int length) {
		MaxLength = e.MaxLength;
		Length = length;
		Charge = e.Charge;

		Mz = new double[MaxLength];
		Intensity = new float[MaxLength];
		memcpy(Mz, e.Mz + start, length * sizeof(double));
		memcpy(Intensity, e.Intensity + start, length * sizeof(float));

		Position = new int[MaxLength];
		memcpy(Position, e.Position + start, length * sizeof(int));

		for (int i = length; i < MaxLength; i++) {
			Mz[i] = 0;
			Intensity[i] = 0;
			Position[i] = -1;
		}
	}


	Rscore* JudgePattern(int element, float normalCutoff = RATIOCUTOFF, float seCutoff = SECUTOFF) {
		Rscore* r = NULL;
		int p = GetMaxIntensityPosition();

		double mass = Mz[GetMaxIntensityPosition()] * Charge;
		Peptide se(mass - ExactMass[element], 1), normal(mass);
		R[0] = Rate(se, Length, ISOLENGTH - Length, element);
		R[1] = Rate(normal, Length, ISOLENGTH - Length, element);

		EnvelopeParameters Parameters{ 0,0,0 };
		Parameters.Length = Length;
		Parameters.Mass = ((int)round(Mz[GetMaxIntensityPosition()] * Charge / 100)) * 100;
		Parameters.Start = R[0].Position;
		if ( R[0].R < seCutoff 
			&& R[1].R / R[0].R > normalCutoff 
			&& AcceptedParameters[element].find(Parameters) != AcceptedParameters[element].end()) {
			r = &R[0];
			Tag &= 0b01111111111111111111111111111111;
		}

		return r;
	}

	Rscore* JudgeEnvelope(float intensity[], int length) {
		R[0] = Rate(intensity, min(static_cast<int>(Length), length), 0);
		return &R[0];
	}

	Rscore* JudgeMeta() {
		double mass = Mz[GetMaxIntensityPosition()] * Charge;
		Peptide meta(mass);
		R[0] = Rate(meta, Length, 0, Se);
		return &R[0];
	}

	bool JudgeMS2Pair(int massShift) {
		int p1 = -1;
		int p2 = -1;

		float* intensitycopy = new float[Length + 2]();
		
		memcpy(intensitycopy + 1, Intensity, Length * sizeof(float));

		for (int i = 1; i <= Length; ++i) {
			if (intensitycopy[i] > intensitycopy[i - 1] && intensitycopy[i] > intensitycopy[i + 1]) {
				if (p1 == -1) {
					p1 = i;
				}
				else {
					p2 = i;
					break;
				}
			}
		}
		if (p2 - p1 <= massShift + 1 
			&& p2 - p1 >= massShift - 1 
			&& min(intensitycopy[p1], intensitycopy[p2]) / max(intensitycopy[p1], intensitycopy[p2]) > 0.7) {
			delete[] intensitycopy;
			return true;
		
		}
		delete[] intensitycopy;
		return false;
	}

	Envelope* Filter(float ratio) {
		double maxIntens = GetMaxIntensity();
		int start = 0;
		for (int i = 0; i < Length; ++i) {
			if (Intensity[i] > maxIntens * ratio) {
				start = i;
				break;
			}
		}

		return new Envelope(*this, start, Length - start);
	}


private:
	unsigned short MaxLength = 0;
	double* Mz = NULL;
	float* Intensity = NULL;
	short Charge;
	unsigned short Length = 0;
	Rscore R[2];

	int* Position = NULL;

	Rscore Rate(Peptide& peptide, int length, int maxShift, int element) {
		Rscore r{ -1 };
		float isodis[ISOLENGTH];
		peptide.CalcDistribution(isodis, element);

		float z = 0;
		for (int j = 0; j < length; j++) {
			z += Intensity[j] * Intensity[j];
		}

		for (int i = 0; i <= maxShift; i++) {
			float x = 0, y = 0;
			//y=0;z=0;
			for (int j = 0; j < length; j++) {
				x += Intensity[j] * isodis[i + j];
				y += isodis[i + j] * isodis[i + j];
			}
			float temp = acos(x / (sqrt(y) * sqrt(z))) * 180.0 / PI;
			if (r.R == -1 || temp < r.R) {
				r.R = temp;
				r.Position = i;
			}
		}

		for (int i = 0; i < Length; i++) {
			r.IsotopeDistribution[i] = isodis[i + r.Position];
		}
		return r;
	}

	Rscore Rate(float isodis[], int length, int maxShift) {
		Rscore r{ -1 };

		float z = 0;
		for (int j = 0; j < length; j++) {
			z += Intensity[j] * Intensity[j];
		}

		for (int i = 0; i <= maxShift; i++) {
			float x = 0, y = 0;
			for (int j = 0; j < length; j++) {
				x += Intensity[j] * isodis[i + j];
				y += isodis[i + j] * isodis[i + j];
			}
			float temp = acos(x / (sqrt(y) * sqrt(z))) * 180.0 / PI;
			if (r.R == -1 || temp < r.R) {
				r.R = temp;
				r.Position = i;
			}
		}

		for (int i = 0; i < Length; i++) {
			r.IsotopeDistribution[i] = isodis[i + r.Position];
		}
		return r;
	}
};

//void FindEnvelope(vector<Envelope*>& envelopes, double mz[], double intensity[], int length, int lenMin = 6, int lenMax = 14) {
//	vector<Envelope*> temp;
//	temp.reserve(length * 1.5f);
//	Envelope* e = NULL;
//	unsigned short charge;
//	for (int k = 0; k < AcceptChargeLength; ++k) {
//		charge = AcceptCharge[k];
//		temp.clear();
//		
//		for (int i = 0; i < length; i++) {
//			if (mz[i] == -1) continue;
//			bool ff = false;
//			double mass = mz[i] * charge;
//			for (int j = temp.size() - 1; j >= 0; j--) {
//				if (temp[j] == NULL) continue;
//				double lastmass = (*temp[j]).GetLastMz() * (*temp[j]).GetCharge();
//				double last2mass = (*temp[j]).GetMz((*temp[j]).GetLength() - 2) * (*temp[j]).GetCharge();
//				if (abs(lastmass + NEUTRON - mass) < ENVELOPEPPM * (mass + lastmass) / 2.0) {
//					(*temp[j]).Add(mz[i], intensity[i], i);
//					temp.push_back(temp[j]);
//					temp[j] = NULL;
//					ff = true;
//				}
//				else if (abs(last2mass + NEUTRON - mass) < ENVELOPEPPM * (mass + last2mass) / 2.0) {
//					e = new Envelope(*temp[j], 0, (*temp[j]).GetLength() - 1);
//					(*e).Add(mz[i], intensity[i], i);
//					temp.push_back(e);
//					ff = true;
//				}
//				else if (lastmass < (mass - NEUTRON) * (1.0 - ENVELOPEPPM)) {
//					break;
//				}
//			}
//			if (!ff) {
//				e = new Envelope(mz[i], intensity[i], i, charge);
//				temp.push_back(e);
//			}
//		}
//		e = NULL;
//		
//		int counttemp = 0;
//		for (int j = temp.size() - 1; j >= 0; j--) {
//			if (temp[j] == NULL) continue;
//			counttemp += 1;
//
//			int len = (*temp[j]).GetLength();
//			bool addflag = false;
//			if (len >= lenMin && len <= lenMax) {
//				addflag = true;
//				for (vector<Envelope*>::iterator k = envelopes.begin(); k != envelopes.end(); k++) {
//					bool containflag = true;
//					for (int jj = 0; jj < len; jj++) {
//						bool uniqueflag = true;
//						for (int kk = 0; kk < (**k).GetLength(); kk++) {
//							if ((*temp[j]).GetMz(jj) == (**k).GetMz(kk)) {
//								uniqueflag = false;
//								break;
//							}
//						}
//						if (uniqueflag) {
//							containflag = false;
//							break;
//						}
//					}
//					if (containflag) {
//						addflag = false;
//						break;
//					}
//				}
//			}
//			if (!addflag) {
//				delete temp[j];
//				temp[j] = NULL;
//			}
//			else {
//				envelopes.push_back(temp[j]);
//			}
//
//		}
//	}
//}
void FindEnvelope(vector<Envelope*>& envelopes, double mz[], double intensity[], int length, int lenMin = 6, int lenMax = 14) {
    vector<unique_ptr<Envelope>> temp;  // Use unique_ptr for automatic cleanup
    temp.reserve(length * 1.5f);
    Envelope* e = nullptr;
    unsigned short charge;
    
    for (int k = 0; k < AcceptChargeLength; ++k) {
        charge = AcceptCharge[k];
        temp.clear();
        
        for (int i = 0; i < length; i++) {
            if (mz[i] == -1) continue;
            bool ff = false;
            double mass = mz[i] * charge;
            
            for (int j = temp.size() - 1; j >= 0; j--) {
                if (!temp[j]) continue;
                
                double lastmass = temp[j]->GetLastMz() * temp[j]->GetCharge();
                double last2mass = temp[j]->GetMz(temp[j]->GetLength() - 2) * temp[j]->GetCharge();
                
                if (abs(lastmass + NEUTRON - mass) < ENVELOPEPPM * (mass + lastmass) / 2.0) {
                    temp[j]->Add(mz[i], intensity[i], i);
                    temp.push_back(move(temp[j]));  // Move ownership
                    ff = true;
                }
                else if (abs(last2mass + NEUTRON - mass) < ENVELOPEPPM * (mass + last2mass) / 2.0) {
                    e = new Envelope(*temp[j], 0, temp[j]->GetLength() - 1);
                    e->Add(mz[i], intensity[i], i);
                    temp.push_back(unique_ptr<Envelope>(e));
                    ff = true;
                }
                else if (lastmass < (mass - NEUTRON) * (1.0 - ENVELOPEPPM)) {
                    break;
                }
            }
            
            if (!ff) {
                temp.push_back(make_unique<Envelope>(mz[i], intensity[i], i, charge));
            }
        }
        
        for (auto& env : temp) {
            if (!env) continue;
            
            int len = env->GetLength();
            bool addflag = (len >= lenMin && len <= lenMax);
            
            if (addflag) {
                for (auto& existing : envelopes) {
                    bool containflag = true;
                    for (int jj = 0; jj < len; jj++) {
                        bool uniqueflag = true;
                        for (int kk = 0; kk < existing->GetLength(); kk++) {
                            if (env->GetMz(jj) == existing->GetMz(kk)) {
                                uniqueflag = false;
                                break;
                            }
                        }
                        if (uniqueflag) {
                            containflag = false;
                            break;
                        }
                    }
                    if (containflag) {
                        addflag = false;
                        break;
                    }
                }
            }
            
            if (addflag) {
                envelopes.push_back(env.release());  // Transfer ownership
            }
        }
    }
}

void FindPatternEnvelope(int element,
                         vector<Envelope*>& seEnvelopes,
                         const vector<Envelope*>& envelopes,
                         float normalCutoff /*= RATIOCUTOFF*/,
                         float seCutoff    /*= SECUTOFF*/,
                         int lenMin        /*= 6*/,
                         int lenMax        /*= 14*/) {
    for (Envelope* env : envelopes) {
        if (!env) {
            continue;
        }

        // 长度范围过滤
        int length = env->GetLength();
        if (length < lenMin || length > lenMax) {
            continue;
        }

        // m/z * 电荷阈值过滤
        double weightedMz = env->GetMaxMz() * env->GetCharge();
        double threshold = 100.0 + ExactMass[element];
        if (weightedMz <= threshold) {
            continue;
        }

        // 调用 JudgePattern，注意它返回的是 Rscore*（或 auto*）
        auto* score = env->JudgePattern(element, normalCutoff, seCutoff);
        if (!score) {
            continue;
        }

        // 找到匹配，收集当前 Envelope*
        seEnvelopes.push_back(env);
    }
}

BYTE* Report(int* count, double mz[], double intensity[], int length, float normalCutoff = RATIOCUTOFF, float seCutoff = SECUTOFF, int lenmin = 6, int lenmax = 14) {
	vector<Envelope*> envelopes;
	cout << "Start FindEnvelope ... ";
	FindEnvelope(envelopes, mz, intensity, length, lenmin, lenmax);
	cout << "Envelope:" << envelopes.size() << endl;

	vector<Envelope*> seEnvelopes;
	cout << "Start FindPattern ... ";
	FindPatternEnvelope(Se, seEnvelopes, envelopes, normalCutoff, seCutoff, lenmin, lenmax);
	cout << "Pattern:" << seEnvelopes.size() << endl;

	*count = seEnvelopes.size();
	BYTE* r = NULL;

	if (*count != 0) {
		r = new BYTE[*count * (REPORTLENGTH + sizeof(bool))];
		BYTE* temp = r;

		for (vector<Envelope*>::iterator i = seEnvelopes.begin(); i != seEnvelopes.end(); i++) {
			memcpy(temp, (**i).Report(), REPORTLENGTH);
			temp += REPORTLENGTH + sizeof(bool);
			if ((**i).GetR(0).R < 10) {
				*(temp - sizeof(bool)) = 1;
			}	
			else {
				*(temp - sizeof(bool)) = 0;
			}
				
		}

		temp = NULL;
	}

	for (vector<Envelope*>::iterator i = envelopes.begin(); i != envelopes.end(); i++) {
		delete* i;
		*i = NULL;
	}

	return r;
}

BYTE* ReportEnvelope(int* count, double mz[], double intensity[], int length, float normalCutoff = RATIOCUTOFF, float seCutoff = SECUTOFF, int lenmin = 3, int lenmax = 14) {
	vector<Envelope*> envelopes;
	FindEnvelope(envelopes, mz, intensity, length, lenmin, lenmax);

	sort(envelopes.begin(), envelopes.end(), [](Envelope* a, Envelope* b) -> bool { return a->GetMaxIntensity() > b->GetMaxIntensity(); });

	*count = envelopes.size();

	BYTE* r = NULL;

	if (*count != 0) {
		r = new BYTE[*count * (REPORTLENGTH + sizeof(bool))];
		BYTE* temp = r;

		for (vector<Envelope*>::iterator i = envelopes.begin(); i != envelopes.end(); i++) {
			memcpy(temp, (**i).Report(), REPORTLENGTH);
			temp += REPORTLENGTH + sizeof(bool);

			*(temp - sizeof(bool)) = 0;
		}

		temp = NULL;
	}

	for (vector<Envelope*>::iterator i = envelopes.begin(); i != envelopes.end(); i++) {
		delete* i;
		*i = NULL;
	}

	return r;
}

void Difference(int length, int mass, int element, float factor = 0.2) {
	Peptide pattern((double)mass - ExactMass[element], 1);
	Peptide normal(mass);

	float isodisPattern[ISOLENGTH];
	pattern.CalcDistribution(isodisPattern, element);

	float isodisNormal[ISOLENGTH];
	normal.CalcDistribution(isodisNormal, element);


	for (; length <= ISOLENGTH; ++length) {
		Rscore r = Rscore{ -1,-1 };
		for (int i = 0; i < ISOLENGTH - length; ++i) {
			double outmax = 0, inmin = 1;
			for (int j = 0; j < ISOLENGTH; ++j) {
				if (j < i || j >= i + length) outmax = max(static_cast<double>(isodisPattern[j]), outmax);
				else inmin = min(static_cast<double>(isodisPattern[j]), inmin);
			}
			if (inmin * (1.0 + factor) < outmax * (1.0 - factor)) continue;


			float z = 0;
			for (int j = 0; j < length; ++j) {
				z += isodisPattern[i + j] * isodisPattern[i + j];
			}

			int diff = ISOLENGTH - length;
			for (int j = 0; j <= diff; ++j) {
				float x = 0, y = 0;
				for (int k = 0; k < length; ++k) {
					x += isodisPattern[i + k] * isodisNormal[j + k];
					y += isodisNormal[j + k] * isodisNormal[j + k];
				}

				float temp = acos(x / (sqrt(y) * sqrt(z))) * 180.0 / PI;
				temp = isnan(temp) ? 0 : temp;
				if (r.R == -1 || temp < r.R) {
					r.R = temp;
					r.Position = j;
				}
			}

			CutoffMap[EnvelopeParameters{ length, i, mass }] = r.R;
			//cout << mass << " " << length << " " << i << " " << r.R << endl;
		}
	}

}

void CalcCutoff(int massMin, int massMax, int lengthMin) {
	massMin = massMin / 100 * 100;
	massMax = ceil(massMax / 100) * 100;

	AcceptedParameters.clear();

	for (int iii = 0; iii < ElementNumber - 5; ++iii) {
		CutoffMap.clear();

		set<EnvelopeParameters> temp;
		for (int mass = massMin; mass <= massMax; mass += 100) {
			Difference(lengthMin, mass, iii);
		}

		for (const auto &p : CutoffMap) {
			if (p.second > InitialCutoff[iii]) {
				temp.insert(p.first);
				//cout << "mass: " << p.first.Mass << " pos:" << p.first.Start << " length:" << p.first.Length << endl;
			}
		}
		AcceptedParameters.push_back(temp);
	}
	
}

int Initialize(int massrangelow = 300, int massrangehigh = 5000) {
	ifstream infile("params.dat", ios::in | ios::binary);
	if (!infile) {
		CalcCutoff(massrangelow, massrangehigh, 5);
		return 1;
	}

	infile.read((char*)&AcceptChargeLength, sizeof(int));
	infile.read((char*)AcceptCharge, sizeof(unsigned short) * AcceptChargeLength);

	for (int i = 0; i < AcceptChargeLength; i++) {
		cout << AcceptCharge[i] << ",";
	}
	cout << endl;

	AcceptedParameters.clear();
	EnvelopeParameters p;
	float elementiso[ELEMENTISOLENGTH];
	float cutoff;
	double mass;
	while (infile.read((char*)& p, sizeof(p))) {
		if (p.Length == -1) {
			infile.read((char*)elementiso, sizeof(float) * ELEMENTISOLENGTH);
			infile.read((char*)&mass, sizeof(double));
			infile.read((char*)&cutoff, sizeof(float));
			AddElement(elementiso, mass, cutoff);
		}
		else {
			(AcceptedParameters.end() - 1)->insert(p);
		}
	}

	//for (auto pp : AcceptedParameters)
	//{
	//	cout << "special:::" << endl;
	//	for (const auto &p : pp) cout << "mass: " << p.Mass << " pos:" << p.Start << " length:" << p.Length << endl;
	//}
		
	infile.close();

	return 1;
}


int AddElementExtern(float distribution[ELEMENTISOLENGTH], double exactmass, float initialcutoff) {
	AddElement(distribution, exactmass, initialcutoff);
	return 1;
}

int ReInitialize(int massrangelow = 300, int massrangehigh = 5000, int minLength = 5) {
	CalcCutoff(massrangelow, massrangehigh, minLength);
	return 1;
}

int SaveParams() {
	ofstream outfile("params.dat", ios::out | ios::binary);
	outfile.write((char*)&AcceptChargeLength, sizeof(int));
	outfile.write((char*)AcceptCharge, sizeof(unsigned short) * AcceptChargeLength);

	EnvelopeParameters startp{ -1,-1,-1 };

	for (int iii = 0; iii < ElementNumber - 5; ++iii) {
		outfile.write((char*)&startp, sizeof(startp));
		outfile.write((char*)(Element + 5 + iii), sizeof(float) * ELEMENTISOLENGTH);
		outfile.write((char*)(ExactMass + iii), sizeof(double));
		outfile.write((char*)(InitialCutoff + iii), sizeof(float));
		for (const auto& p : AcceptedParameters[iii]) {
			outfile.write((char*)&p, sizeof(p));
		}
	}
	outfile.close();
	return 1;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <ms1_data.txt>\n";
        return EXIT_FAILURE;
    }

    ifstream infile(argv[1]);
    if (!infile) {
        cerr << "Error: cannot open file " << argv[1] << "\n";
        return EXIT_FAILURE;
    }

    //Init
    Initialize();

    vector<double> mzs;
    vector<double> ints;
    string line;
    int spectrum_idx = 0;

    auto process_spectrum = [&]() {
        if (mzs.empty()) return;
        ++spectrum_idx;

        cout << "Spectrum " << spectrum_idx
                  << ": points=" << mzs.size()
                  << endl;

        int count = 0;
        unsigned char* results = Report(&count, mzs.data(), ints.data(), static_cast<int>(mzs.size()));
        cout << "count:" << count << endl;

        delete [] results;
        mzs.clear();
        ints.clear();
    };

    while (getline(infile, line)) {
        if (line.empty()) {
            continue;
        }
        // 如果是新谱开始标志，先处理上一谱
        if (line[0] == 'S') {
            process_spectrum();
            continue;
        }
        // 跳过所有以 'I' 开头的注释行
        if (line[0] == 'I') {
            continue;
        }
        // 否则应为数据行：四列数值，取前两列 mz 和 intensity
        istringstream iss(line);
        double mz, intensity;
        // 读两列，忽略后面两列
        if (iss >> mz >> intensity) {
            mzs.push_back(mz);
            ints.push_back(intensity);
        }
    }
    // 文件末尾，处理最后一谱
    process_spectrum();

    delete[] AcceptCharge;
    delete[] Element;
    delete[] ExactMass;
    delete[] InitialCutoff;

    return EXIT_SUCCESS;
}

