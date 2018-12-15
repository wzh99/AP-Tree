#include <iostream>
#include <ctime>
#include "../include/aptree.hpp";

using namespace std;
Query generateQuery(const string k1, const string k2);

int main() 
{
	vector<Query> queries;

	srand(time(NULL));
	queries.push_back(generateQuery("Bruceville", "Cemetery	"));
	queries.push_back(generateQuery("Boulder", "Agua"));
	queries.push_back(generateQuery("Sal", "Creek"));
	queries.push_back(generateQuery("Agua", "Sal"));
	queries.push_back(generateQuery("Wash", "Aguaje"));
	queries.push_back(generateQuery("Draw", "Valley"));
	queries.push_back(generateQuery("Arlington", "State"));
	queries.push_back(generateQuery("Wildlife", "Area"));
	queries.push_back(generateQuery("Park", "Bar"));

	vector<string> vocab;
	vocab.push_back("Agua");
	vocab.push_back("Aguaje");
	vocab.push_back("Area");
	vocab.push_back("Arlington");
	vocab.push_back("Bar");
	vocab.push_back("Boulder");
	vocab.push_back("Bruceville");
	vocab.push_back("Cemetery");
	vocab.push_back("Creek");
	vocab.push_back("Draw");
	vocab.push_back("Park");
	vocab.push_back("Sal");
	vocab.push_back("State");
	vocab.push_back("Valley");
	vocab.push_back("Wash");
	vocab.push_back("Wildlife");

	STObject object1, object2;
	object1 = genereteObject("Sal", "Agua", "Wildlife");
	object2 = genereteObject("Bourder", "Draw", "Park");
}

Query generateQuery(const string k1, const string k2){
	Query obj;
	obj.keywords.insert(k1);
	obj.keywords.insert(k2);

	double x1, x2, y1, y2, x11, x22, y11, y22;
	x1 = (double)rand() / (double)RAND_MAX;
	x2 = (double)rand() / (double)RAND_MAX;
	y1 = (double)rand() / (double)RAND_MAX;
	y2 = (double)rand() / (double)RAND_MAX;
	obj.region.min.x = (x1 < x2) ? x1 : x2;//This may not be able to generate uniform distribution.
	obj.region.min.y = (y1 < y2) ? y1 : y2;//But I think uniform distribution is not necessary so I keep it :)
	obj.region.max.x = (x1 < x2) ? x2 : x1;
	obj.region.max.y = (y1 < y2) ? y2 : y1;
	
	return obj;
}

STObject genereteObject(const string s1, const string s2, const string s3) {
	STObject obj;
	obj.keywords.insert(s1);
	obj.keywords.insert(s2);
	obj.keywords.insert(s3);

	double xLocation = (double)rand() / (double)RAND_MAX;
	double yLocation = (double)rand() / (double)RAND_MAX;
	obj.location.x = xLocation;
	obj.location.y = yLocation;

	return obj;
}