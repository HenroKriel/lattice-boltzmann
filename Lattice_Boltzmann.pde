//https://physics.weber.edu/schroeder/javacourse/LatticeBoltzmann.pdf //<>// //<>//
//left click to disturb the fluid, right click to add obstacles!
//the screen is colored according to the speed of the fluid

import java.util.*;

final int grid_len = 200;

//the discrete velocity vectors (directions that the "particles" can move)
final PVector[] e = new PVector[9];

//the weights corresponding to the directions
final float[] w = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};

//controls viscosity of fluid
float visc = 1;

//number densities of particles moving in the 9 different directions.
float[][][] n = new float[9][grid_len][grid_len];

//total density at each lattice point
float[][] p = new float[grid_len][grid_len];

//average velocity at each lattice point
PVector[][] u = new PVector[grid_len][grid_len];

//physical objects that the particles bounce off
boolean[][] obstacles = new boolean[grid_len][grid_len];

//higher means more sensitive
float sensitivity = 15;

void setup() {
  size(800, 800, FX2D);

  e[0] = new PVector(0, 0);
  e[1] = new PVector(1, 0);
  e[2] = new PVector(0, 1);
  e[3] = new PVector(-1, 0);
  e[4] = new PVector(0, -1);
  e[5] = new PVector(1, 1);
  e[6] = new PVector(-1, 1);
  e[7] = new PVector(-1, -1);
  e[8] = new PVector(1, -1);

  //creates boarder of obstacles
  for (int i = 0; i < grid_len; i++) {
    obstacles[i][0] = true;
    obstacles[i][grid_len-1] = true;
    obstacles[0][i] = true;
    obstacles[grid_len-1][i] = true;
  }

  //fill the room with fluid
  for (int j = 1; j < grid_len - 1; j++) {
    for (int k = 1; k < grid_len - 1; k++) {
      n[0][j][k] = 1;
    }
  }
}

void draw() {
  background(255);

  //puts right facing fountain in middle of room
  n[1][grid_len/2][grid_len/2] = 1;

  if (mousePressed) {
    float grid_size = width/grid_len;
    int x = (int) (mouseX/grid_size);
    int y = (int) (mouseY/grid_size);
    
    //creates disturbance
    if (mouseButton == LEFT) {
      for (int i = 0; i < 9; i++) {
        n[i][x][y] = 1;
      }
    }
    
    if (mouseButton == RIGHT) {
      for(int i = -1; i <= 1; i++) {
        for(int j = -1; j <= 1; j++) {
          obstacles[x+i][y+j] = true;
        }
      }
    }
  }

  collisions();
  streaming();

  draw_particles();
  draw_obstacles();

  //draw_grid();
}

void collisions() {
  //clears p and u to be recalculated
  for (int j = 0; j < grid_len; j++) {
    for (int k = 0; k < grid_len; k++) {
      p[j][k] = 0;
      u[j][k] = new PVector(0, 0);
    }
  }

  //calculates p and u
  for (int j = 0; j < grid_len; j++) {
    for (int k = 0; k < grid_len; k++) {
      for (int i = 0; i < 9; i++) {
        p[j][k] += n[i][j][k];
        u[j][k].add(PVector.mult(e[i], n[i][j][k]));
      }
      if (p[j][k] != 0) u[j][k].div(p[j][k]);
      
      //prevents model from becoming unstable
      //(so far has only happened when room is not filled with fluid)
      u[j][k].limit(sqrt(6)/3);
    }
  }

  //this is the collision step
  for (int j = 0; j < grid_len; j++) {
    for (int k = 0; k < grid_len; k++) {
      for (int i = 0; i < 9; i++) {
        float n_eq = ( p[j][k]*w[i] * (1 + 3*PVector.dot(e[i], u[j][k]) + 9.0/2*pow(PVector.dot(e[i], u[j][k]), 2) - 3.0/2*pow(u[j][k].mag(), 2)) );
        float n_old = n[i][j][k];
        n[i][j][k] = n_old + visc*(n_eq-n_old);
      }
    }
  }
}

void streaming() {
  //we need this separate memory so that moving the particles doesn't affect later calculations
  float[][][] new_n = new float[9][grid_len][grid_len];

  //moves particles to lattice point in respective direction
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < grid_len; j++) {
      for (int k = 0; k < grid_len; k++) {
        int x2 = ( j + (int) e[i].x );
        int y2 = ( k + (int) e[i].y );

        //bounds checking
        if (x2 >= 0 && x2 < grid_len && y2 >= 0 && y2 < grid_len) {
          //particles are "bounced back" when they hit an obstacle, reversing their direction
          if (obstacles[x2][y2]) {
            int reverse_i;
            //I want to cycle between directions (1, 2, 3, 4) or (5, 6, 7, 8)
            //I figured out that (i-shift + step)%4 + shift gives the desired result
            //step is 2 (because the opposite direction of 1 is 3, etc)
            //shift is the start of the cycle, i.e. 1, 5
            if (i == 0) reverse_i = 0;
            else if (i <= 4) reverse_i = (i+1)%4 + 1;
            else reverse_i = (i-3)%4 + 5;
            
            new_n[reverse_i][j][k] += n[i][j][k];
          } else {
            new_n[i][x2][y2] += n[i][j][k];
          }
        }
      }
    }
  }

  n = new_n;
}

void color_particle(int i, int j, float x) {

  //sets x to value between 0 and 1

  //caps x off after max_particles
  //final float max_particles = 0.05;
  //x = constrain(x, 0, max_particles)/max_particles;

  //sigmoid function gives a lot resolution to smaller values
  //without removing variation in higher values like method above
  x = sensitivity*x/sqrt(1+pow(sensitivity*x, 2));

  color background = color(120, 0, 0);
  color fluid = color(255, 100, 100);
  color c = lerpColor(background, fluid, x);
  fill_square(i, j, c);
}

void fill_square(int i, int j, color c) {
  push();
  noStroke();
  fill(c);
  rect(i*(width/grid_len), j*(height/grid_len), width/grid_len, height/grid_len);
  pop();
}

void draw_grid() {
  for (int i = 0; i < grid_len; i++) {
    line(i*(width/grid_len), 0, i*(width/grid_len), height);
    line(0, i*(height/grid_len), width, i*(height/grid_len));
  }
}

void draw_particles() {
  for (int j = 0; j < grid_len; j++) {
    for (int k = 0; k < grid_len; k++) {
      //colors according to speed. Can be changed to density (p[j][k])
      color_particle(j, k, u[j][k].mag());
    }
  }
}

void draw_obstacles() {
  for (int j = 0; j < grid_len; j++) {
    for (int k = 0; k < grid_len; k++) {
      if (obstacles[j][k])
        fill_square(j, k, color(0, 0, 100));
    }
  }
}
