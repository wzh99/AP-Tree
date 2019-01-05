#include <iostream>
#include <fstream>
#include <set>
#include <ctime>
#include "../include/aptree.hpp"
#include "clock.hpp"

using namespace std;
void generateTestData(int querySize, int queryNumber, int objectSize, int objectNumber, int vocabSize,
					  double maxQueryEdge, double minQueryEdge,
					  vector<Query> &query, vector<STObject> &object, vector<string> &vocab);
void varyDataParameters(int querySize, int queryNumber, double maxQueryEdge, double minQueryEdge);
void varyConstructParameters(int f, int ��q);

//test���ݣ� number of queries 1000~5000,����ʱ�䣬�Լ�����object��matchʱ��ı仯
			//query�����Ӱ�죻
			//query�ؼ�������Ӱ�죻

int main()
{
	//effect of f
	cout << "f:\n";
	//varyConstructParameters(16, 10);
	//for (int varyF = 10; varyF <= 100; varyF+=10)
	//	varyConstructParameters(varyF, 10);

	//effect of ��q
	cout << "��q:\n";
	for (int vary��q = 10; vary��q < 300; vary��q += 10) varyConstructParameters(16, vary��q);
}

void varyConstructParameters(int f, int ��q) {
	vector<Query> query;
	vector<STObject> object;
	vector<string> vocab;

	generateTestData(3, 2000, 7, 2000, 1000, 0.1, 0.01, query, object, vocab);

	APTree *tree;
	//cout << "Construction parameters: f:" << f << " ��q:" << ��q << "\nconstruction: ";
//	TIME_COUNT(
		tree = new APTree(vocab, query, f, ��q)
		;
//			, 1, milliseconds);

	//cout << "match 2000: ";
	TIME_COUNT(
		for (int i = 0; i < 2000; i++) tree->Match(object[i])
			, 1, microseconds);

	delete tree;
}

void varyDataParameters(int querySize, int queryNumber, double maxQueryEdge, double minQueryEdge) {
	vector<Query> query;
	vector<STObject> object;
	vector<string> vocab;

	generateTestData(querySize, queryNumber, 7, 20000, 10000, maxQueryEdge, minQueryEdge, query, object, vocab);

	APTree *tree;
	cout << "Data parameters: querySize:" << querySize << " queryNumber:" << queryNumber 
		<< " maxQueryEdge:" << maxQueryEdge << " minQueryEdge:" << minQueryEdge << "\nconstruction: ";
	TIME_COUNT(tree = new APTree(vocab, query, 16, 10), 1, milliseconds);

	cout << "match 20000: ";
	TIME_COUNT(
		for (int i = 0; i < 20000; i++) tree->Match(object[i])
	, 1, microseconds);

	delete tree;
}

void generateTestData(int querySize, int queryNumber, int objectSize, int objectNumber, int vocabSize, 
					  double maxQueryEdge, double minQueryEdge, 
					  vector<Query> &query,  vector<STObject> &object, vector<string> &vocab){

	if (vocabSize > 10000) { cout << "\nVocab Size is too large!"; return; }
	vocab.clear();
	object.clear();
	query.clear();

	ifstream testDataFile("NationalFileWordOrdered.txt");
	string wordTmp;
	vector<string> wholeVocabulary;
	for (int i = 0; i < 10286; i++) {
		testDataFile >> wordTmp;
		wholeVocabulary.push_back(wordTmp);
	}
	
	int indexTmp;
	for (int VS = 0; VS < vocabSize; VS++) {
		indexTmp = rand() % wholeVocabulary.size();
		vocab.push_back(wholeVocabulary[indexTmp]);
		wholeVocabulary.erase(wholeVocabulary.begin() + indexTmp);
	}
	
	for (int i = 0; i < objectNumber; i++) {
		STObject objectTmp;
		set<string> tmpObjectSet;

		while (tmpObjectSet.size() < objectSize)
			tmpObjectSet.insert(vocab[rand() % vocabSize]);//put randomly chosen keywords in order and get rid of identical keywords

		while (tmpObjectSet.size() > 0) {
			objectTmp.keywords.insert(*tmpObjectSet.begin());
			tmpObjectSet.erase(tmpObjectSet.begin());
		}

		objectTmp.location.x = (double)rand() / (double)RAND_MAX;
		objectTmp.location.y = (double)rand() / (double)RAND_MAX;

		object.push_back(objectTmp);
	}
	
	for (int i = 0; i < queryNumber; i++) {
		int ActualQuerySize = (querySize > 0) ? querySize : (rand() % 5 + 1);
		//if user gives negative or zero value to queryNumber, then the size of every query will be randomly chosen between 1 and 5
		Query queryTmp;
		set<string> tmpQuerySet;
		while(tmpQuerySet.size() < ActualQuerySize)
			tmpQuerySet.insert(vocab[rand() % vocabSize]);//put randomly chosen keywords in order and get rid of identical keywords
		while(tmpQuerySet.size() > 0) {
			queryTmp.keywords.insert(*tmpQuerySet.begin());
			tmpQuerySet.erase(tmpQuerySet.begin());
		}
		double xTraverse = ((double)rand() / (double)RAND_MAX) * (maxQueryEdge - minQueryEdge) + minQueryEdge;
		double yTraverse = ((double)rand() / (double)RAND_MAX) * (maxQueryEdge - minQueryEdge) + minQueryEdge;

		queryTmp.region.min.x  = ((double)rand() / (double)RAND_MAX) * (1 - xTraverse);
		queryTmp.region.min.y  = ((double)rand() / (double)RAND_MAX) * (1 - yTraverse);
		queryTmp.region.max.x = queryTmp.region.min.x + xTraverse;
		queryTmp.region.max.y = queryTmp.region.min.y + yTraverse;

		query.push_back(queryTmp);
	}
}

//Query generateQuery(const string k1, const string k2);
//STObject generateObject(const string s1, const string s2, const string s3);
//
//int main() 
//{
//	vector<Query> queries;
//
//	srand(time(NULL));
//	queries.push_back(generateQuery("Bruceville", "Cemetery"));
//	queries.push_back(generateQuery("Boulder", "Agua"));
//	queries.push_back(generateQuery("Sal", "Creek"));
//	queries.push_back(generateQuery("Aguaje", "Sal"));
//	queries.push_back(generateQuery("Wash", "Aguaje"));
//	queries.push_back(generateQuery("Draw", "Valley"));
//	queries.push_back(generateQuery("Arlington", "State"));
//	queries.push_back(generateQuery("Wildlife", "Area"));
//	queries.push_back(generateQuery("Park", "Bar"));
//
//	vector<string> vocab;
//	vocab.push_back("Agua");
//	vocab.push_back("Aguaje");
//	vocab.push_back("Area");
//	vocab.push_back("Arlington");
//	vocab.push_back("Bar");
//	vocab.push_back("Boulder");
//	vocab.push_back("Bruceville");
//	vocab.push_back("Cemetery");
//	vocab.push_back("Creek");
//	vocab.push_back("Draw");
//	vocab.push_back("Park");
//	vocab.push_back("Sal");
//	vocab.push_back("State");
//	vocab.push_back("Valley");
//	vocab.push_back("Wash");
//	vocab.push_back("Wildlife");
//	cout << generateQuery("Bruceville", "Cemetery");
//
//	STObject object1, object2;
//	object1 = generateObject("Sal", "Creek", "Wildlife");
//	object2 = generateObject("Boulder", "Draw", "Park");
//	cout << object1;
//
//	APTree aptree(vocab, queries, 4, 4);
//	aptree.Match(object1);
//	aptree.Match(object2);
//	return 0;
//}
//
//Query generateQuery(const string k1, const string k2){
//	Query obj;
//	obj.keywords.insert(k1);
//	obj.keywords.insert(k2);
//
//	double x1, x2, y1, y2;
//	x1 = (double)rand() / (double)RAND_MAX;
//	x2 = (double)rand() / (double)RAND_MAX;
//	y1 = (double)rand() / (double)RAND_MAX;
//	y2 = (double)rand() / (double)RAND_MAX;
//	obj.region.min.x = (x1 < x2) ? x1 : x2;
//	obj.region.min.y = (y1 < y2) ? y1 : y2;
//	obj.region.max.x = (x1 < x2) ? x2 : x1;
//	obj.region.max.y = (y1 < y2) ? y2 : y1;
//	
//	return obj;
//}
//
//STObject generateObject(const string s1, const string s2, const string s3) {
//	STObject obj;
//	obj.keywords.insert(s1);
//	obj.keywords.insert(s2);
//	obj.keywords.insert(s3);
//
//	double xLocation = (double)rand() / (double)RAND_MAX;
//	double yLocation = (double)rand() / (double)RAND_MAX;
//	obj.location.x = xLocation;
//	obj.location.y = yLocation;
//
//	return obj;
//}