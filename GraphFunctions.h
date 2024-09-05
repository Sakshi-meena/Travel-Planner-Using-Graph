#ifndef GRAPHFUNCTIONS_H
#define GRAPHFUNCTIONS_H

#include <string>
#include <vector>
#include <iostream>
#include <queue>
#include <cmath>      // Include cmath for fabs
#include <stack>

#include "Location.h"
#include "Route.h"
#include "FileOperations.h"

// Define compareLocation struct if not defined elsewhere
struct compareLocation {
    bool operator()(Location* l1, Location* l2) {
        return l1->lengthFromStart > l2->lengthFromStart;  // Min-heap based on lengthFromStart
    }
};

class Graph {
public:
    std::vector<Location*> cities;
    std::vector<Route*> routes;

    int numExists;

    Graph(const std::string& nodesFile, const std::string& edgesFile);

    int getCityIndex(const std::string& key) const;
    Location* getCity(const std::string& key) const;

    float getWeight(const std::string& startS, const std::string& endS, bool costOrTime) const;
    float getWeight(Location* start, Location* end, bool costOrTime) const;

    void Dijkstras(const std::string& startS, bool costOrTime);

    std::vector<Location*>* adjacentLocations(Location* city) const;
    Route* getRoute(Location* start, bool costOrTime, float totalDistance) const;

    std::vector<Route*>* adjacentRoutes(Location* city) const;
    
    std::stack<Location*> cityStacker(const std::string& destinationS) const;
    std::stack<Route*> routeStacker(const std::string& destinationS, bool costOrTime) const;

    //-----------------------------------------------------------------------
    /*	SCRAPPED METHODS
    
        Location* getSmallest();
        void printShortestRoute(const std::string& destinationS);
        void printOutRoutes();
    */
    //-----------------------------------------------------------------------
};

Graph::Graph(const std::string& nodesFile, const std::string& edgesFile) {
    routes = routeParser(edgesFile);
    cities = locationParser(nodesFile, routes);

    numExists = cities.size();
}

int Graph::getCityIndex(const std::string& key) const {
    for (int i = 0; i < cities.size(); ++i) {
        if (cities[i]->country == key) {
            return i;
        }
    }
    return -1;
}

Location* Graph::getCity(const std::string& key) const {
    int i = getCityIndex(key);
    return (i != -1) ? cities[i] : nullptr;
}

float Graph::getWeight(const std::string& startS, const std::string& endS, bool costOrTime) const {
    Location* start = getCity(startS);
    Location* end = getCity(endS);

    for (const auto& route : routes) {
        if (route->doesConnect(start, end)) {
            return costOrTime ? route->cost : route->time;
        }
    }

    return -1;
}

float Graph::getWeight(Location* start, Location* end, bool costOrTime) const {
    for (const auto& route : routes) {
        if (route->doesConnect(start, end)) {
            return costOrTime ? route->cost : route->time;
        }
    }

    return -1;
}

void Graph::Dijkstras(const std::string& startS, bool costOrTime) {
    Location* start = getCity(startS);
    start->lengthFromStart = 0;

    std::priority_queue<Location*, std::vector<Location*>, compareLocation> minHeap;

    for (const auto& city : cities) {
        minHeap.push(city);
    }

    while (!minHeap.empty()) {
        while (!minHeap.empty() && !minHeap.top()->exists) {
            minHeap.pop();
        }

        if (minHeap.empty()) {
            return;
        }

        Location* smallest = minHeap.top();
        minHeap.pop();

        smallest->exists = false;

        auto adjacentCities = adjacentLocations(smallest);

        for (const auto& adjacent : *adjacentCities) {
            float distance = getWeight(smallest, adjacent, costOrTime) + smallest->lengthFromStart;

            if (distance < adjacent->lengthFromStart) {
                adjacent->lengthFromStart = distance;
                adjacent->previous = smallest;
            }

            minHeap.push(adjacent);  // Reinsert into heap with updated distance
        }

        delete adjacentCities;
    }
}

std::vector<Location*>* Graph::adjacentLocations(Location* city) const {
    auto output = new std::vector<Location*>();

    for (const auto& route : city->routes) {
        if (route->destination->exists) {
            output->push_back(route->destination);
        }
    }

    return output;
}

std::vector<Route*>* Graph::adjacentRoutes(Location* city) const {
    auto output = new std::vector<Route*>();

    for (const auto& route : routes) {
        if (route->origin->capital == city->capital) {
            output->push_back(route);
        }
    }

    return output;
}

Route* Graph::getRoute(Location* start, bool costOrTime, float totalDistance) const {
    auto adjacentRoutes = this->adjacentRoutes(start);

    const float epsilon = 1e-5;

    for (const auto& route : *adjacentRoutes) {
        float routeMeasure = costOrTime ? route->cost : route->time;
        if (std::fabs((totalDistance - routeMeasure) - route->origin->lengthFromStart) > epsilon) {
            delete adjacentRoutes;
            return route;
        }
    }

    delete adjacentRoutes;
    return nullptr;
}

std::stack<Location*> Graph::cityStacker(const std::string& destinationS) const {
    Location* destination = getCity(destinationS);
    std::stack<Location*> output;

    while (destination) {
        output.push(destination);
        destination = destination->previous;
    }

    return output;
}

std::stack<Route*> Graph::routeStacker(const std::string& destinationS, bool costOrTime) const {
    std::stack<Route*> output;
    Location* destination = getCity(destinationS);
    float totalDistance = destination->lengthFromStart;

    while (destination->previous) {
        output.push(getRoute(destination->previous, costOrTime, totalDistance));
        destination = destination->previous;
        totalDistance = destination->lengthFromStart;
    }

    return output;
}

/*
Location* Graph::getSmallest() {
    if (cities.empty()) {
        return nullptr;
    }

    Location* output = cities[0];
    for (int i = 1; i < cities.size(); ++i) {
        if (cities[i]->lengthFromStart < output->lengthFromStart) {
            output = cities[i];
        }
    }

    cities.erase(std::remove(cities.begin(), cities.end(), output), cities.end());
    --numExists;
    return output;
}

void Graph::printShortestRoute(const std::string& destinationS) {
    Location* destination = getCity(destinationS);
    Location* previous = destination;
    float distFromStart = 0;
    while (previous) {
        std::cout << previous->capital << " <- ";
        distFromStart += previous->lengthFromStart;
        previous = previous->previous;
    }
    std::cout << std::endl;
    std::cout << "Distance from start: " << distFromStart << std::endl;
}

void Graph::printOutRoutes() {
    for (const auto& city : cities) {
        std::cout << "Routes from: " << city->capital << std::endl;
        for (const auto& route : city->routes) {
            std::cout << "    " << route->destination->capital << " Type: " << route->transport << " Cost: " << route->cost << " Time: " << route->time << std::endl;
        }
    }
}
*/

#endif
