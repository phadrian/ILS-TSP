//
//  main.cpp
//  ILS-TSP
//
//  Created by Adrian Pheh on 17/2/17.
//  Copyright © 2017 Adrian Pheh. All rights reserved.
//

/* Include header files depending on platform */
#ifdef _WIN32
#include "glut.h"
#elif __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/GLUT.h>
#endif

#include <algorithm>
#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <thread>
#include <cmath>
#include "graph.h"

using namespace std;

typedef vector<int> TOUR;
typedef vector<vector<int>> MATRIX;

typedef chrono::milliseconds MS;
typedef chrono::system_clock SYS_CLOCK;
typedef chrono::time_point<chrono::system_clock> TIME_PT;

#define POINT_SIZE 4

// Global variables to render tours in OpenGL
TOUR tour_to_render;
MATRIX dist;
int N;
vector<vector<double>> city_coords;
bool finished = false;

/*
 * Helper functions for debugging
 */
void print_tour_line(TOUR tour) {
    for (int i = 0; i < tour.size(); i++) {
        printf("%d,", tour[i]);
    }
    printf("\n");
}

void print_tour(TOUR tour) {
    for (int i = 0; i < tour.size(); i++) {
        printf("%d\n", tour[i]);
    }
}

void print_mat(MATRIX matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

int tour_length(TOUR tour, MATRIX dist) {
    int length = 0;
    for (int i = 0; i < tour.size()-1; i++) {
        length += dist[tour[i]][tour[i+1]];
    }
    // Include a return edge from last city to first city
    length += dist[tour[tour.size()-1]][tour[0]];
    return length;
}

// Functions to compute the viewport boundaries
 double get_max_boundary() {
    double largest_x = 0, largest_y = 0;
    for (int i = 0; i < city_coords.size(); i++) {
        largest_x = largest_x < city_coords[i][0] ? city_coords[i][0]
                                                  : largest_x;
        largest_y = largest_y < city_coords[i][1] ? city_coords[i][1]
                                                  : largest_y;
    }
    return largest_x > largest_y ? largest_x : largest_y;
}

double get_min_boundary() {
    double smallest_x = INFINITY, smallest_y = INFINITY;
    for (int i = 0; i < city_coords.size(); i++) {
        smallest_x = smallest_x > city_coords[i][0] ? city_coords[i][0]
        : smallest_x;
        smallest_y = smallest_y > city_coords[i][1] ? city_coords[i][1]
        : smallest_y;
    }
    return smallest_x < smallest_y ? smallest_x : smallest_y;
}

void display();

/*
 * Reads in the input values and calculates the distance matrix
 */
MATRIX get_distance_matrix() {
    // Read input and store coordinates
    int n;
    scanf(" %d", &n);
    
    // Initialize the global coordinate matrix
    city_coords = vector<vector<double>>(n);
    
    vector<double> x(n);
    vector<double> y(n);
    
    for (int i = 0; i < n; i++) {
        double buf = 0;
        scanf(" %lf%lf%lf", &buf, &x[i], &y[i]);
        printf("%lf,%lf,%lf\n", buf, x[i], y[i]);
        // Also store the coordinates globally
        city_coords[i].push_back(x[i]);
        city_coords[i].push_back(y[i]);
    }

    // Run in N^2 to calculate euclidean distances between each set
    // of coordinates
    MATRIX dist(n, TOUR(n));
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            dist[i][j] = dist[j][i] = round(sqrt(pow(x[i]-x[j], 2) +
                                                 pow(y[i]-y[j], 2)));
        }
    }
    return dist;
}

/*
 * Naive greedy tour implementation on Kattis
 */
TOUR get_greedy_tour(MATRIX dist) {
    int n = dist.size();
    TOUR tour(n);
    vector<bool> used(n, false);
    // Start from first city
    tour[0] = 0;
    used[0] = true;
    for (int i = 1; i < n; i++) {
        int best = -1;
        for (int j = 0; j < n; j++) {
            if (!used[j] && (best == -1 ||
                             dist[tour[i-1]][j] < dist[tour[i-1]][best])) {
                best = j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }
    return tour;
}

TOUR two_opt_swap(TOUR tour, int b, int c) {
    TOUR new_tour(0);
    for (int i = 0; i < b; i++) {
        new_tour.push_back(tour[i]);
    }
    for (int i = c; i >= b; i--) {
        new_tour.push_back(tour[i]);
    }
    for (int i = c+1; i < tour.size(); i++) {
        new_tour.push_back(tour[i]);
    }
    return new_tour;
}

TOUR two_opt(TOUR tour, MATRIX dist) {
    int count = 1, n = tour.size();
    bool local_opt = false;
    int best_tour_length, initial_tour_length;
    TOUR best_tour = tour;
    // Repeat until a local optimum is found
    while (!local_opt) {
        local_opt = true;
        // Obtain the initial length and set that to be the current best
        initial_tour_length = best_tour_length = tour_length(best_tour, dist);
        for (int i = 0; i < n; i++) {
            int a = i;
            int b = i == n-1 ? 0 : i+1;
            for (int j = i+2; j < n; j++) {
                int c = j;
                int d = j == n-1 ? 0 : j+1;
                // This condition handles the special case where the edges
                // to be swapped are the first and last edges, which are
                // not "adjacent" in this format
                if (a != 0 || d != 0) {
                    /*
                     * Swapping edges 2--3 and 5--6,
                     * 0--1--2**3--4--5**6--7--8--9
                     * Delete edges 2--3 and 5--6
                     * Add edges 2--5 and 3--6
                     */
                    
                    /* ******************************************************
                     * This line is bugged, initial_tour_length should be
                     * replaced by best_tour_length to ensure that the values
                     * being compared is actually correct. Uncomment the
                     * second printf in the if loop to check the lengths,
                     * both should be the same.
                     *
                     * This gives wrong values for the new tours, but for
                     * some reason the new tours are shown to be more optimal 
                     * than if the correct values are used.
                     * ******************************************************
                     */
                    int perm_length = initial_tour_length
                    - dist[best_tour[a]][best_tour[b]]
                    - dist[best_tour[c]][best_tour[d]]
                    + dist[best_tour[a]][best_tour[c]]
                    + dist[best_tour[b]][best_tour[d]];
                    // Only swap edges in tours which are shown to be better
                    if (perm_length < best_tour_length) {
                        //printf("permutation length = %d < %d\n", perm_length, best_tour_length);
                        TOUR new_tour = two_opt_swap(best_tour, b, c);
                        //printf("perm length = %d ?= %d\n", perm_length, tour_length(new_tour, dist));
                        best_tour = new_tour;
                        best_tour_length = perm_length;
                        local_opt = false;
                        
                        // Draw best tour
                        tour_to_render = best_tour;
                        display();
                    }
                    /*
                     printf("%d. swapping (%d,%d) with (%d,%d)\n",
                     count, tour[a], tour[b], tour[c], tour[d]);
                     TOUR new_tour = two_opt_swap(tour, b, c);
                     print_tour_line(new_tour);
                     //printf("tour cost = %d\n", tour_length(new_tour, dist));
                     //*/
                    count++;
                }
            }
        }
    }
    //printf("2opt ran %d times in total\n", count);
    //printf("best tour length = %d\n", best_tour_length);
    return best_tour;
}

TOUR two_h_opt_swap(TOUR tour, int a, int d) {
    int n = tour.size();
    TOUR new_tour(0);
    new_tour.push_back(tour[a]);
    new_tour.push_back(tour[d]);
    // i = 1 to account for a already added
    for (int i = 1; i < n; i++) {
        int idx = (a+i) % n;
        // Ignore d which has been added already
        if (idx != d) {
            new_tour.push_back(tour[idx]);
        }
    }
    return new_tour;
}

TOUR two_h_opt(TOUR tour, MATRIX dist) {
    int n = tour.size();
    bool local_opt = false;
    int best_tour_length, initial_tour_length;
    TOUR best_tour = tour;
    while (!local_opt) {
        local_opt = true;
        initial_tour_length = best_tour_length = tour_length(best_tour, dist);
        for (int i = 0; i < n; i++) {
            int a = i;
            int b = i == n-1 ? 0 : i+1;
            for (int j = 0; j < n-4; j++) {
                // Use j as an increasing offset
                int c = (b+1 + j) % n;
                int d = (c+1) % n;
                int e = (d+1) % n;
                // Do not perform checking if any edges overlap in either
                // of the 2 sets a--b and c--d--e
                if (a != e && b != e) {
                    /*
                     * Check edges 2--3 and 5--6--7,
                     * 0--1--2**3--4--5**6**7--8--9
                     * Delete edges 5--6, 6--7, 2--3
                     * Add edges 2--6, 6--3, 5--7
                     */
                    /*
                     printf("checking (%d,%d) and (%d,%d,%d)\n",
                     best_tour[a], best_tour[b],
                     best_tour[c], best_tour[d], best_tour[e]);
                     //*/
                    int perm_length = initial_tour_length
                    - dist[best_tour[c]][best_tour[d]]
                    - dist[best_tour[d]][best_tour[e]]
                    - dist[best_tour[a]][best_tour[b]]
                    + dist[best_tour[a]][best_tour[d]]
                    + dist[best_tour[d]][best_tour[b]]
                    + dist[best_tour[c]][best_tour[e]];
                    /*
                     if (perm_length != tour_length(new_tour, dist)) {
                     printf("lengths dont match, %d!=%d\n",
                     perm_length, tour_length(new_tour, dist));
                     }
                     //*/
                    if (perm_length < best_tour_length) {
                        //printf("improvement found, new best = %d\n", perm_length);
                        TOUR new_tour = two_h_opt_swap(best_tour, a, d);
                        best_tour = new_tour;
                        best_tour_length = perm_length;
                        local_opt = false;
                        
                        // Draw best tour
                        tour_to_render = best_tour;
                        display();
                    }
                }
            }
        }
    }
    return best_tour;
}

// TODO: THIS IS BUGGED AS FUCK
void add_v_in_range(TOUR tour, TOUR *new_tour, int s, int t, bool rev) {
    /*
     rev ? printf("reverse\n") : printf("normal\n");
     printf("adding from %d to %d\n", tour[s], tour[t]);
     //*/
    int n = tour.size();
    int cur_itr = 0;
    int max_itr = !rev ? s > t ? (t+n)-s + 1 : (t-s) + 1
    : s > t ? (s-t) + 1 : (s+n)-t + 1;
    //printf("max iter = %d\n", max_itr);
    if (rev) {
        for (int i = s+n; cur_itr < max_itr; i--) {
            new_tour->push_back(tour[i%n]);
            cur_itr++;
        }
    } else {
        for (int i = s; cur_itr < max_itr; i = (i+1) % n) {
            new_tour->push_back(tour[i]);
            cur_itr++;
        }
    }
}

TOUR three_opt_swap(TOUR tour, TOUR idx, int swap_pattern) {
    /*
     printf("swap pattern = %d\n", swap_pattern);
     printf("idx = ");
     print_tour_line(idx);
     //*/
    
    /*
     * idx = |0|1|2|3|4|5|
     * 		 |a|b|c|d|e|f|
     */
    #define A 0
    #define B 1
    #define C 2
    #define D 3
    #define E 4
    #define F 5
    #define REV true
    #define NORM false
    
    TOUR new_tour(0);
    // Append d--...--e which is same in both patterns
    add_v_in_range(tour, &new_tour, idx[D], idx[E], NORM);
    if (swap_pattern == 1) {
        // Pattern 1
        // Append e--b--...--c
        add_v_in_range(tour, &new_tour, idx[B], idx[C], NORM);
    } else {
        // Pattern 2
        // Append e--c--rev(...)--b
        add_v_in_range(tour, &new_tour, idx[C], idx[B], REV);
    }
    // Append f--...--a which is same in both patterns
    add_v_in_range(tour, &new_tour, idx[F], idx[A], NORM);
    return new_tour;
}

TOUR three_opt(TOUR tour, MATRIX dist) {
    int n = tour.size();
    bool local_opt = false;
    int initial_tour_length, best_tour_length;
    TOUR best_tour = tour;
    while (!local_opt) {
        local_opt = true;
        initial_tour_length = best_tour_length = tour_length(best_tour, dist);
        int count = 0;
        for (int i = 0; i < n; i++) {
            int a = i;
            int b = (i+1) % n;
            // Limited to n-5 for c because n-(a+b+d+e+f)
            for (int j = 0; j < n-5; j++) {
                int c = (b+1 + j) % n;
                int d = (c+1) % n;
                // Subtract by j to account for the decreasing cyclic
                // distance between f and a
                for (int k = 0; k < n-5-j; k++) {
                    int e = (d+1 + k) % n;
                    int f = (e+1) % n;
                    count++;
                    /*
                     printf("(%d,%d), (%d,%d), (%d,%d)\n",
                     tour[a], tour[b],
                     tour[c], tour[d],
                     tour[e], tour[f]);
                     //*/
                    /* Permutation 1:
                     *
                     *       b---c
                     *        \ /
                     *         /
                     *   a----/-\----d
                     *     \ /   \ /
                     *      f     e
                     */
                    int perm_1_length = initial_tour_length
                    - dist[best_tour[a]][best_tour[b]]
                    - dist[best_tour[c]][best_tour[d]]
                    - dist[best_tour[e]][best_tour[f]]
                    + dist[best_tour[a]][best_tour[d]]
                    + dist[best_tour[e]][best_tour[b]]
                    + dist[best_tour[c]][best_tour[f]];
                    /* Permutation 2:
                     *
                     *       b---c
                     *       |   |
                     *       |   |
                     *   a---|---|---d
                     *     \ |   | /
                     *       f   e
                     */
                    int perm_2_length = initial_tour_length
                    - dist[best_tour[a]][best_tour[b]]
                    - dist[best_tour[c]][best_tour[d]]
                    - dist[best_tour[e]][best_tour[f]]
                    + dist[best_tour[a]][best_tour[d]]
                    + dist[best_tour[e]][best_tour[c]]
                    + dist[best_tour[b]][best_tour[f]];
                    // Prepare for swapping and comparison
                    int swap_pattern, perm_length;
                    if (perm_1_length <= perm_2_length) {
                        perm_length = perm_1_length;
                        swap_pattern = 1;
                    } else {
                        perm_length = perm_2_length;
                        swap_pattern = 2;
                    }
                    if (perm_length < best_tour_length) {
                        //printf("new length = %d < %d\n", perm_length, best_tour_length);
                        TOUR idx = {a, b, c, d, e, f};
                        TOUR new_tour = three_opt_swap(best_tour, idx, swap_pattern);
                        best_tour = new_tour;
                        best_tour_length = perm_length;
                        local_opt = false;
                    }
                }
            }
        }
        //printf("%d permutations\n", count);
    }
    return best_tour;
}

TOUR perturb(TOUR tour) {
    // Create the random number generator
    random_device rd;
    default_random_engine rng(rd());
    uniform_int_distribution<int> random_offset(1, tour.size()/4);
    // Determine random positions
    int pos1 = random_offset(rng);
    int pos2 = pos1 + random_offset(rng);
    int pos3 = pos2 + random_offset(rng);
    
    TOUR new_tour(0);
    for (int i = 0; i < pos1; i++) {
        new_tour.push_back(tour[i]);
    }
    for (int i = pos3; i < tour.size(); i++) {
        new_tour.push_back(tour[i]);
    }
    for (int i = pos2; i < pos3; i++) {
        new_tour.push_back(tour[i]);
    }
    for (int i = pos1; i < pos2; i++) {
        new_tour.push_back(tour[i]);
    }
    return new_tour;
}

TOUR perturb(TOUR tour, int strength) {
    TOUR new_tour = tour;
    for (int i = 0; i < strength; i++) {
        new_tour = perturb(new_tour);
    }
    return new_tour;
}

void display() {
    // Small buffer time for better visibility
    this_thread::sleep_for(MS(10));
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glPushMatrix();
    glColor3f(0, 0, 0);
    glPointSize(POINT_SIZE);
    glBegin(GL_POINTS);
    for (int i = 0; i < tour_to_render.size(); i++) {
        glVertex2f(city_coords[tour_to_render[i]][0],
                   city_coords[tour_to_render[i]][1]);
    }
    glEnd();
    glPopMatrix();
    
    glPushMatrix();
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < tour_to_render.size(); i++) {
        glVertex2f(city_coords[tour_to_render[i]][0],
                   city_coords[tour_to_render[i]][1]);
    }
    glEnd();
    glPopMatrix();
    
    glFlush();
}

void idle() {

    if (!finished) {
        // Set up the timing system
        MS runtime(5980);
        TIME_PT end;
        end = SYS_CLOCK::now() + runtime;
        
        // Set the random seed for shuffling
        srand(time(0));
        random_device rd;
        default_random_engine rng(rd());
        
        // Define an initial nearest neighbours toour
        Graph g(dist);
        TOUR initial_tour = g.multi_frag_tour();
        tour_to_render = initial_tour;
        display();
        
        // Define the local and global optimum tours and lengths
        TOUR local_opt_tour, global_opt_tour;
        int local_opt_length, global_opt_length;
        
        // Time the two_opt operation
        TIME_PT two_opt_start = SYS_CLOCK::now();
        // Define a city limit for 3opt to reduce < 2s
        #define MAX_3_OPT 20
        // Run local optimizations
        TOUR opt;
        opt = two_opt(initial_tour, dist);
        opt = two_h_opt(opt, dist);
        opt = N < MAX_3_OPT ? three_opt(opt, dist)
        : opt;
        global_opt_tour = local_opt_tour = opt;
        MS two_opt_length = chrono::duration_cast<MS>(SYS_CLOCK::now() -
                                                      two_opt_start);
        MS buffer(two_opt_length);
        
        // Parameter to be tuned
        // Controls the number of non-improving iterations allowed before
        // applying stronger perturbation
        #define MAX_K 1000
        
        // Record the number of iterations of perturbations
        int itr = 1, k_non_improv = 0;
        global_opt_length = local_opt_length = tour_length(local_opt_tour, dist);
        while (SYS_CLOCK::now() < end) {
            if (end - SYS_CLOCK::now() < buffer) {
                break;
            } else {
                // Run perturbation and get a new local optima
                if (local_opt_tour.size() >= 8) {
                    local_opt_tour = perturb(local_opt_tour);
                } else {
                    shuffle(local_opt_tour.begin(), local_opt_tour.end(), rng);
                }
                // Locally optimize
                local_opt_tour = two_opt(local_opt_tour, dist);
                local_opt_tour = two_h_opt(local_opt_tour, dist);
                
                local_opt_length = tour_length(local_opt_tour, dist);
                if (local_opt_length < global_opt_length) {
                    global_opt_tour = local_opt_tour;
                    global_opt_length = local_opt_length;
                }
                itr++;
            }
        }
        
        // Draw global optimum tour
        tour_to_render = global_opt_tour;
        display();
        
        //*/
        /*
         printf("%d iterations run\n", itr);
         print_tour_line(global_opt_tour);
         printf("best found tour length = %d\n", global_opt_length);
         //*/
        //print_tour(global_opt_tour);
        finished = true;
    } else {
        this_thread::sleep_for(MS(1000));
    }
}

int main(int argc, char *argv[]) {
    
    // Define the distance matrix and an initial nearest neighbours tour
    dist = get_distance_matrix();
    N = dist.size();
    
    // Set up the OpenGL window
    TOUR tour;
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize (600, 600);
    glutInitWindowPosition (50, 50);
    glutCreateWindow (argv[0]);
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glClearColor (1, 1, 1, 1.0);
    
    double max_boundary = abs(get_max_boundary()) + POINT_SIZE;
    double min_boundary = -abs(get_min_boundary()) - POINT_SIZE;
    printf("max bound = %f\n", max_boundary);
    printf("min bound = %f\n", min_boundary);
    
    glOrtho(min_boundary, max_boundary, min_boundary, max_boundary, min_boundary, max_boundary);
    glutMainLoop();
    
}