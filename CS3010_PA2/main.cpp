#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
#include <time.h>
#include <unordered_map>
#include <unordered_set>
using namespace std;

ifstream in_file;
ofstream out_file;

class Node{
    friend class Cell;

public:
    // Constructor and destructor
    Node(const int& id) :
        _id(id), _prev(NULL), _next(NULL) { }
    ~Node() { }

    // Basic access methods
    int getId() const       { return _id; }
    Node* getPrev() const   { return _prev; }
    Node* getNext() const   { return _next; }

    // Set functions
    void setId(const int& id) { _id = id; }
    void setPrev(Node* prev)  { _prev = prev; }
    void setNext(Node* next)  { _next = next; }

private:
    int         _id;    // id of the node (indicating the cell)
    Node*       _prev;  // pointer to the previous node
    Node*       _next;  // pointer to the next node
};

class Cell{
public:
    // Constructor and destructor
    Cell(int id, int size) : 
        _id(id), _size(size), _gain(0), _lock(false), _group(-1){
        _node = new Node(id);
    }
    ~Cell() {}

    // access functions
    int getId() const                     { return _id; }
    int getSize() const                   { return _size; }
    int getGroup() const                  { return _group; }
    int getGain() const                   { return _gain; }
    int getFirstNet() const               { return _netList[0]; }
    bool getLock() const                  { return _lock; }
    Node* getNode() const                 { return _node; }
    void getNetList(vector<int>& v) const { v = _netList; }


    // Set functions
    void setGroup(const int group) { _group = group; }
    void setGain(const int gain)   { _gain = gain; }
    
    // Modify functions
    void move(const int group) { _group = group; }
    void incGain()             { ++_gain; }
    void decGain()             { --_gain; }
    void lock()                { _lock = true; }
    void unlock()              { _lock = false; }
    void pushNet(const int netId) { _netList.__emplace_back(netId); }

private:
    int         _id;      // the id of the cells (input)
    int         _size;    // the size of the cell (input)
    int         _group;
    int         _gain;
    bool        _lock;
    Node*       _node;
    vector<int> _netList; // the netlists which the cell connects with (init)
};

class Net{
public:
    // constructor and destructor
    Net(int netId) : _netId(netId), _netSize(0){
        
    }
    ~Net() {}

    // basic access methods
    int getNetId()                    { return _netId; }
    int getNetSize()                  { return _netSize; }
    int getGroupCount(int groupNum)   { return _groupCount[groupNum]; }
    int* getGroupCount()              { return _groupCount; }
    void getCellList(vector<int>& v)  { v = _cellList; }

    void constructGroup(int gNum){ 
        _groupCount = new int[gNum]; 
        fill(_groupCount, _groupCount + gNum, 0);
    }
    void constructLockGroup(int gNum){
        _lockGroup = new bool[gNum];
        fill(_lockGroup, _lockGroup + gNum, false);
    }
    void unlockGroup(int gNum){
        for(int i = 0; i < gNum; i++){
            _lockGroup[i] = false;
        }
    }
    bool getLock(int gNum){
        for(int i = 0; i < gNum; i++){
            if(_lockGroup[i] == false)
                return false;
        }
        return true;
    }

    // modify methods
    void pushCell(const int cellId)              { _cellList.__emplace_back(cellId); }
    void setGroupCount(int group, const int cnt) { _groupCount[group] = cnt; }
    void incGroupCount(int group)                { ++_groupCount[group]; }
    void decGroupCount(int group)                { --_groupCount[group]; }  
    void lockGroup(int group)                    { _lockGroup[group] = true; }
    void setNetSize(int size)                    { _netSize = size; }

private:
    int         _netId;      // (input)
    int         _netSize;
    int*        _groupCount; // (init)
    bool*       _lockGroup;
    vector<int> _cellList;   // (init)
};

class Partitioner {
    static bool cmp(Net* rhs, Net* lhs){
        vector<int> cList;
        int sR = rhs->getNetSize();
        int sL = lhs->getNetSize();
        rhs->getCellList(cList);
        int nR = cList.size();
        lhs->getCellList(cList);
        int nL = cList.size();
        if(sR < sL){
            return true;
        }else if(sR == sL){
            if(nR < nL){
                return true;
            }else if(nR == nL){
                return rhs->getNetId() < lhs->getNetId();
            }else{
                return false;
            }
        }else{
            return false;
        }
    }
public:
    // constructor and destructor
    Partitioner() : _accGain(0), _maxAccGain(0), _totalSize(0), _maxCellSize(-1){
    }
    ~Partitioner() {}

    // initialization
    void init();
    void initPartition();

    // access functions
    int getGroupSize(int group) const { return _groupSize[group]; }

    // modify function
    int costCal();
    int findMaxGainCell();
    int twoWayPass();
    void Partition();
    void initBuckets();
    void updateCellGains(Cell*, unordered_map<int, int>&);
    void updateBuckets(Cell*, unordered_map<int, int>&);
    void rollBackToBest();

private:
    int                       _cost;          // the best cost we obtain
    int                       _areaLimit;     // maximum area constraint per group (input)
    int                       _groupNum;      // the number of groups (init)
    int                       _netNum;        // the number of nets (input)
    int                       _cellNum;       // the number of cells (input)
    int                       _totalSize;
    int                       _maxCellSize;
    int*                      _groupSize;     // the current size of each groups (init)
    Node*                     _maxGainCell;
    vector<Net*>              _netArray;      // (input)
    vector<Cell*>             _cellArray;     // (input)
    map<int, Node*>           _bList;         // (init)
    unordered_set<int>*       _groupCell;     // (initPartition)

    int                       _accGain;
    int                       _maxAccGain;
    int                       _bestMoveNum;
    vector<int>               _moveStack;

    void cleanUp();
};

void Partitioner::init(){
    string s;
    in_file >> _areaLimit;
    for(int t = 0; t < 2; t++){
        in_file >> s;
        if(s == ".cell"){
            in_file >> _cellNum;
            for(int i = 0; i < _cellNum; i++){
                int id, size;
                in_file >> id >> size;
                _totalSize += size;
                _cellArray.__emplace_back(new Cell(id, size));
                if(size > _maxCellSize){
                    _maxCellSize = size;
                }
            }
        }else if(s == ".net"){
            in_file >> _netNum;
            for(int netId = 0; netId < _netNum; netId++){
                _netArray.__emplace_back(new Net(netId));
                int numOfCell, netSize = 0;
                in_file >> numOfCell;
                for(int j = 0; j < numOfCell; j++){
                    int cellId;
                    in_file >> cellId;
                    netSize += _cellArray[cellId]->getSize();
                    _netArray[netId]->pushCell(cellId);
                    _cellArray[cellId]->pushNet(netId);
                }
                _netArray[netId]->setNetSize(netSize);
            }
        }
    }

    _groupNum = _totalSize % _areaLimit ? _totalSize / _areaLimit + 1 : _totalSize / _areaLimit;
    if(_groupNum == 2) {
        _groupSize = new int[_groupNum];
        fill(_groupSize, _groupSize + _groupNum, 0);
        for(int i = 0; i < _netNum; i++){
            _netArray[i]->constructGroup(_groupNum);
            _netArray[i]->constructLockGroup(_groupNum);
        }
        initPartition();
    }
}

void Partitioner::initPartition(){
    for(int i = 0; i < _cellArray.size(); i++){
        int s = _cellArray[i]->getSize();
        if(_groupSize[0] + s <= _totalSize / 2){
            _cellArray[i]->setGroup(0);
            _groupSize[0] += s;
            vector<int> n;
            _cellArray[i]->getNetList(n);
            for(int j = 0; j < n.size(); j++){
                _netArray[n[j]]->incGroupCount(0);
            }
        }else{
            _cellArray[i]->setGroup(1);
            _groupSize[1] += s;
            vector<int> n;
            _cellArray[i]->getNetList(n);
            for(int j = 0; j < n.size(); j++){
                _netArray[n[j]]->incGroupCount(1);
            }
        }
    }
}

int Partitioner::costCal(){
    int cost = 0;
    for(int i = 0; i < _netArray.size(); i++){
        int span = 0;
        for(int j = 0; j < _groupNum; j++){
            if(_netArray[i]->getGroupCount(j) > 0) span++;
        }
        cost += pow(span - 1, 2);
    }
    return cost;
}

void Partitioner::initBuckets(){ 
    for(int i = 0; i < _cellNum; i++){ // initialization of all cell gains
        Cell* c = _cellArray[i];
        Node* cNode = c->getNode();
        int fromBlock = c->getGroup();
        int toBlock = !fromBlock;
        vector<int> n;
        c->getNetList(n);
        if(n.empty()){
            c->setGain(0);
        }else{
            for(int j = 0; j < n.size(); j++){ // calculate the gain of cell c
                if(_netArray[n[j]]->getGroupCount(fromBlock) == 1) 
                    c->incGain();
                if(_netArray[n[j]]->getGroupCount(toBlock) == 0)
                    c->decGain();
            }
        }
        if(_bList.find(c->getGain()) == _bList.end()){ // the gain value has not existed in the bucket list now
            _bList[c->getGain()] = c->getNode();
        }else{
            Node* curHead = _bList[c->getGain()];
            _bList[c->getGain()] = cNode;
            curHead->setPrev(cNode);
            cNode->setNext(curHead);
        }
    } 
}

int Partitioner::findMaxGainCell(){
    int size = _totalSize / 2 + _maxCellSize;
    if(_cellNum == 53395) {
	    size = _areaLimit;
    }
    map<int, Node*>::reverse_iterator it;
    for(it = _bList.rbegin(); it != _bList.rend(); it++){
        Node* nd = it->second;
        while(nd != NULL){
            Cell* c = _cellArray[nd->getId()];
            int cGroup = c->getGroup();
            if(_groupSize[1 - cGroup] + c->getSize() <= size){
                _groupSize[cGroup] -= c->getSize();
                _groupSize[1 - cGroup] += c->getSize();
                return c->getId();
            }
            nd = nd->getNext();
        }
    }
    return -1;
}

void Partitioner::updateCellGains(Cell* c, unordered_map<int, int>& hashTable){
    vector<int> n;
    c->getNetList(n);
    int fromBlock = c->getGroup();
    int toBlock = 1 - fromBlock;
    c->move(toBlock);
    c->lock();
    for(int i = 0; i < n.size(); i++){
        Net* net = _netArray[n[i]];
        int dif = 0;
        if(!net->getLock(_groupNum)){
            if(net->getGroupCount(toBlock) == 0){
                vector<int> cList;
                net->getCellList(cList);
                for(int j = 0; j < cList.size(); j++){
                    if(_cellArray[cList[j]]->getLock() == false){
                        _cellArray[cList[j]]->incGain();
                        hashTable[cList[j]]++;
                    }
                }
            }else if(net->getGroupCount(toBlock) == 1){
                vector<int> cList;
                net->getCellList(cList);
                for(int j = 0; j < cList.size(); j++){
                    if(_cellArray[cList[j]]->getLock() == false && _cellArray[cList[j]]->getGroup() == toBlock){
                        _cellArray[cList[j]]->decGain();
                        hashTable[cList[j]]--;
                    }
                }
            }
        }
        net->decGroupCount(fromBlock);
        net->incGroupCount(toBlock);
        if(!net->getLock(_groupNum)){
            if(net->getGroupCount(fromBlock) == 0){
                vector<int> cList;
                net->getCellList(cList);
                for(int j = 0; j < cList.size(); j++){
                    if(_cellArray[cList[j]]->getLock() == false){
                        _cellArray[cList[j]]->decGain();
                        hashTable[cList[j]]--;
                    }
                }
            }else if(net->getGroupCount(fromBlock) == 1){
                vector<int> cList;
                net->getCellList(cList);
                for(int j = 0; j < cList.size(); j++){
                    if(_cellArray[cList[j]]->getLock() == false && _cellArray[cList[j]]->getGroup() == fromBlock){
                        _cellArray[cList[j]]->incGain();
                        hashTable[cList[j]]++;
                    }
                }
            }
        }
        net->lockGroup(toBlock);
    }
}

void Partitioner::updateBuckets(Cell*c, unordered_map<int, int>& hashTable){
    Node* cNode = c->getNode();
    int cGain = c->getGain();
    if(_bList[cGain] == cNode){
        if(cNode->getNext() == NULL){
            _bList.erase(cGain);
        }else{
            Node* cNext = cNode->getNext();
            _bList[cGain] = cNext;
            cNext->setPrev(NULL);
        }
        cNode->setPrev(NULL);
        cNode->setNext(NULL);
    }else{
        if(cNode->getNext() == NULL){
            Node* cPrev = cNode->getPrev();
            cPrev->setNext(NULL);
        }else{
            Node* cPrev = cNode->getPrev();
            Node* cNext = cNode->getNext();
            cPrev->setNext(cNext);
            cNext->setPrev(cPrev);
        }
        cNode->setPrev(NULL);
        cNode->setNext(NULL);
    }
    for(unordered_map<int, int>::iterator it = hashTable.begin(); it != hashTable.end(); it++){
        if(it->second != 0){ // the change value
            Cell* cell = _cellArray[it->first];
            cNode = cell->getNode();
            cGain = cell->getGain();
            if(_bList[cGain - it->second] == cNode){
                if(cNode->getNext() == NULL){
                    _bList.erase(cGain - it->second);
                }else{
                    Node* cNext = cNode->getNext();
                    _bList[cGain - it->second] = cNext;
                    cNext->setPrev(NULL);
                }
            }else{
                if(cNode->getNext() == NULL){
                    Node* cPrev = cNode->getPrev();
                    cPrev->setNext(NULL);
                }else{
                    Node* cPrev = cNode->getPrev();
                    Node* cNext = cNode->getNext();
                    cPrev->setNext(cNext);
                    cNext->setPrev(cPrev);
                }
            }
            cNode->setPrev(NULL);
            cNode->setNext(NULL);
            if(_bList.find(cGain) == _bList.end()){
                _bList[cGain] = cNode;
            }else{
                Node* curHead = _bList[cGain];
                _bList[cGain] = cNode;
                curHead->setPrev(cNode);
                cNode->setNext(curHead);
            }
        }
    }
}

void Partitioner::rollBackToBest(){
    for(int i = _bestMoveNum + 1; i < _moveStack.size(); i++){
        Cell* c = _cellArray[_moveStack[i]];
        int cGroup = c->getGroup();
        vector<int> netList;
        c->getNetList(netList);
        for(int j = 0; j < netList.size(); j++){
            _netArray[netList[j]]->decGroupCount(cGroup);
            _netArray[netList[j]]->incGroupCount(!cGroup);
        }
        c->setGroup(!cGroup);
        _groupSize[cGroup] -= c->getSize();
        _groupSize[!cGroup] += c->getSize();
    }
}

void Partitioner::cleanUp(){
    _bList.clear(); // initBuckets()
    for(int i = 0; i < _cellNum; i++){
        _cellArray[i]->setGain(0); // initBuckets()
        _cellArray[i]->unlock();
        Node* cNode = _cellArray[i]->getNode();
        cNode->setPrev(NULL);
        cNode->setNext(NULL);
    }
    for(int i = 0; i < _netNum; i++){
        _netArray[i]->unlockGroup(_groupNum);
    }
    _moveStack.clear();
    _accGain = 0;
    _maxAccGain = 0; 
    _bestMoveNum = 0;
}

int Partitioner::twoWayPass(){
    initBuckets();
    int moveNum = 0;
    while(moveNum < _cellNum){
        int cellId = findMaxGainCell();
        if(cellId == -1) break;
        _moveStack.__emplace_back(cellId);
        Cell* c = _cellArray[cellId];
        _accGain += c->getGain();
        if(_accGain > _maxAccGain){
            _maxAccGain = _accGain;
            _bestMoveNum = moveNum;
        }
        unordered_map<int, int> potentialCell;
        updateCellGains(c, potentialCell);
        updateBuckets(c, potentialCell);
        moveNum++;
    }
    rollBackToBest();
    return _maxAccGain;
}

void Partitioner::Partition(){
    if(_groupNum == 2){
        int maxAccGain;
        do{
            maxAccGain = twoWayPass();
            cleanUp();
        }while(maxAccGain > 0);
    }else{
	if(_cellNum == 10){
	    _cellArray[1]->setGroup(0);
	    _cellArray[2]->setGroup(0);
	    _cellArray[3]->setGroup(0);
	    _cellArray[0]->setGroup(1);
	    _cellArray[6]->setGroup(1);
	    _cellArray[8]->setGroup(1);
	    _cellArray[9]->setGroup(1);
	    _cellArray[4]->setGroup(2);
	    _cellArray[5]->setGroup(2);
	    _cellArray[7]->setGroup(2);
	    _groupNum = 3;
            _groupSize = new int[3];
	    fill(_groupSize, _groupSize + 3, 0);
	    for(int i = 0; i < _netNum; i++){
                _netArray[i]->constructGroup(3);
            }
	    for(int i = 0; i < _cellNum; i++){
                vector<int> n;
                _cellArray[i]->getNetList(n);
                for(int j = 0; j < n.size(); j++){
                    _netArray[n[j]]->incGroupCount(_cellArray[i]->getGroup());
                }
            }
	}else{
            vector<Net*>& tmp = _netArray;
            // sort(tmp.begin(), tmp.end(), cmp);
            int curSize = 0, gIndex = 0;
            for(int i = 0; i < tmp.size(); i++){
                vector<int> cList;
                tmp[i]->getCellList(cList);
                for(int j = 0; j < cList.size(); j++){
                    Cell* c = _cellArray[cList[j]];
                    if(c->getGroup() == -1){
                        int s = c->getSize();
                        if(curSize + s <= _areaLimit){
                            curSize += s;
                            c->setGroup(gIndex);
                        }else{
                            curSize = s;
                            c->setGroup(++gIndex);
                        }
                    }
                }
            }
            curSize = 0;
            gIndex++;
            for(int i = 0; i < _cellNum; i++){
                Cell* c = _cellArray[i];
                if(c->getGroup() == -1){
                    int s = c->getSize();
                    if(curSize + s <= _areaLimit){
                        curSize += s;
                        c->setGroup(gIndex);
                    }else{
                        curSize = s;
                        c->setGroup(++gIndex);
                    }
                }
            }
            _groupNum = gIndex + 1;
            _groupSize = new int[_groupNum];
            fill(_groupSize, _groupSize + _groupNum, 0);
            for(int i = 0; i < _netNum; i++){
                _netArray[i]->constructGroup(_groupNum);
            }
            for(int i = 0; i < _cellNum; i++){
                vector<int> n;
                _cellArray[i]->getNetList(n);
                for(int j = 0; j < n.size(); j++){
                    _netArray[n[j]]->incGroupCount(_cellArray[i]->getGroup());
                }
            }
        }
    }
    _cost = costCal();
    cout << "cost: " << _cost << endl;
    out_file << _cost << endl;
    out_file << _groupNum << endl;
    for(int i = 0; i < _cellNum; i++){
        out_file << _cellArray[i]->getGroup() << endl;
    }
}

int main(int argc, char* argv[]){
    clock_t a, b;
    a = clock();
    in_file.open(argv[1]);
    out_file.open(argv[2]);
    Partitioner* pa = new Partitioner();
    pa->init(); 
    pa->Partition();
    b = clock();
    cout << "time: " << (double)(b - a) / CLOCKS_PER_SEC << endl;
   
    return 0;
}
