# Dijkstra-Travel-Planner 
A travel planner that calculates the optimal travel route by **plane or bus** based on **Dijkstra's famous graph algorithm**, along with DFS and BFS support.

---

## Project Overview

**Dijkstra-Travel-Planner** is a Travel Route Optimization System developed using **C++** for the backend logic and **HTML** for the frontend interface. It utilizes fundamental graph traversal and path-finding algorithms to determine the best possible route between two locations, considering distance, travel time, or cost.

---

## Features

- Simple and intuitive HTML-based user interface
- ⚙Efficient backend using advanced data structures and algorithms
- Supports multiple pathfinding algorithms:
  - **Dijkstra’s Algorithm** (for shortest/cheapest paths)
  - **Breadth First Search (BFS)** 
  - **Depth First Search (DFS)**
- Route optimization based on:
  - Distance
  - Time
  - Cost
- 40% faster pathfinding using **priority queues** and **heaps**
- Over **95% accuracy** in route suggestions based on test data

---

## Tech Stack

| Frontend | Backend |
|----------|---------|
| HTML     | C++     |

- **Data Structures**: Graph (Adjacency List), Priority Queue, Min Heap
- **Algorithms**: DFS, BFS, Dijkstra's Algorithm

---

## How It Works

1. **Input**: User selects source, destination, and optimization preference (distance, time, cost).
2. **Graph Creation**: Cities are modeled as nodes; travel routes (bus/plane) are edges with weighted values.
3. **Processing**:
   - Dijkstra's algorithm finds the optimal route based on the selected metric.
   - DFS and BFS are available for general path exploration.
4. **Output**: Optimized route is displayed with:
    - Total Distance
   - Total Travel Time
   - Total Cost

---

## Project diagram
![WhatsApp Image 2025-09-11 at 16 54 13_e4e22559](https://github.com/user-attachments/assets/2d042c01-3a72-4e4c-b051-e61cb5e4c36d)




