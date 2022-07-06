double setpt[nact]={0};
int step_no = 0;
int prev_step = 0;
double theta0 = 0;
double theta_ref0 = 0;

double Kp_ik = 2;

int flag_constant = 1; //1 for constant yaw and 0 for varying yaw

double px = 0.05;
double py = 0;
#define t_curve 60
double tend = t_curve;
double a = 2;

#include <sys/time.h>

#include "utility_matrix.c"
#include "set_model_state.c"
#include "poly_trajectory.c"
#include "get_Jac.c"
#include "get_kinematics2.c"
#include "inverse_kinematics_3D_analytic.c"
#include "control_params.h"
#include "get_trq_stance.c"
#include "get_trq_swing.c"
#include "set_limit.c"
#include "leg_state_machine.c"
#include "get_stance_forces.c"
#include "get_trajectory.c"
#include "nlopt.h"


FILE *fidmat,*fidC; //,*fidC2;
char pathmat[] = "../sample/a1_trot19traj/data.m";
char pathC[] = "../sample/a1_trot19traj/data.txt";
int flag_openmat;
int flag_C;
// ************ Change this as per the problem being solved *********//
#define DATA_PTS 3 //Set this based on columns in data file
float data[5000][DATA_PTS]; //Structure that will store the data.
int STEPS;
int COUNTER=0;
// **************************************************************** //

// void read_data()
// {
//     int i, j;

// 	//fid = fopen(DATA_FILE,"r");
//   fidC = fopen(pathC,"r");

// 	/* Read file */
//     i = 0;
//     while(!feof(fidC))
//     {
//         for (j=0;j<DATA_PTS;j++)
//             fscanf(fidC, "%f", &data[i][j]);
//         i = i+1;
//     }
//     fclose(fidC);
//     STEPS = i-1;

// 	/* For Display only */
//     // for (i=0;i<STEPS;i++)
//     // {
//     //     for (j=0;j<DATA_PTS;j++)
//     //         printf("%f ",data[i][j]);
//     //     printf("\n");
//     // }
// }



int fsm[4] = {0};

void init_controller(int motiontime, double quat_act[4], double omega_act[3],
                     double q_act[nact], double u_act[nact])
{
    int i,j;

    sdinit(); sdprinterr(stderr);


    for (j=0;j<4;j++)
      fsm[j] = fsm_init;

    fidmat = fopen(pathmat,"w");
    fprintf(fidmat,"all_data = [ \n");
    flag_openmat = 1;

  fidC = fopen(pathC,"w");
  flag_C = 1;

  double ang1,ang2, theta;
  double dc[3][3];
  sdquat2dc(quat_act[1],quat_act[2],quat_act[3],quat_act[0],dc); sdprinterr(stderr);
  sddc2ang(dc,&ang1,&ang2,&theta); sdprinterr(stderr);
  theta0 = theta;
  theta_ref0 = theta;

}



















// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
/*

  - the most difficult part is loop inside mycost
    controls[vx, vy, omega] -> [model_diffcar] -> state_update[xdot, ydot, theta_dot] -> [i+1]

*/



double mycost(unsigned n, const double *U) //unsigned n, const double *U
{
  // manual iteration initialization --- the total iteration time decided by the receding horizon
  // double XX[3] = {X[0], X[1], X[2]}; // using local variables to avoid error
  double ts = 0.002;
  double XX[3] = {x_p_global, y_p_global, theta0_global};
  // printf("states inside mycost is %f %f %f\n", XX[0],XX[1],XX[2]);
  double Xref[3] = {x_ref_global, y_ref_global, theta0_ref_global};
  // printf("ref inside mycost is %f %f %f\n", Xref[0],Xref[1],Xref[2]);
  double X_Xref[3];                  // X-Xref
  double mul[3];                     // transpose(X-Xref)*Q
  double cost;                       // cost initialization
  double cost_state;                 // transpose(X-Xref)*Q*(X-Xref)
  double A[9];                       // this diff car model doesn't have A matrix
  double B[9];                       // B matrix of the state space model of diff car
  double B_U[3];                     // B*U
  double Xdot[3];                    // state velocity
  double Q[9] = {1,0,0, 0,1,0, 0,0,1};    // quadrutic function parameters can be tuned
  // double U[15] = {};                      // we have the special index format here (counts by column)

  // iteration 01
    X_Xref[0] = XX[0]-Xref[0];
    X_Xref[1] = XX[1]-Xref[1];
    X_Xref[2] = XX[2]-Xref[2];
    mul[0] = X_Xref[0]*Q[0]+X_Xref[1]*Q[3]+X_Xref[2]*Q[6];
    mul[1] = X_Xref[0]*Q[1]+X_Xref[1]*Q[4]+X_Xref[2]*Q[7];
    mul[2] = X_Xref[0]*Q[2]+X_Xref[1]*Q[5]+X_Xref[2]*Q[8];
    cost_state = X_Xref[0]*mul[0]+mul[1]*X_Xref[1]+mul[2]*X_Xref[2];
    cost = cost + cost_state;
    B[0] = cos(XX[2]);
    B[1] = -sin(XX[2]);
    B[2] = 0;
    B[3] = sin(XX[2]);
    B[4] = cos(XX[2]);
    B[5] = 0;
    B[6] = 0;
    B[7] = 0;
    B[8] = 1;
    B_U[0] = B[0]*U[0]+B[1]*U[1]+B[2]*U[2];
    B_U[1] = B[3]*U[0]+B[4]*U[1]+B[5]*U[2];
    B_U[2] = B[6]*U[0]+B[7]*U[1]+B[8]*U[2];
    Xdot[0] = B_U[0];
    Xdot[1] = B_U[1];
    Xdot[2] = B_U[2];
    XX[0] = XX[0]+ts*Xdot[0];
    XX[1] = XX[1]+ts*Xdot[1];
    XX[2] = XX[2]+ts*Xdot[2];

    // iteration 02
    X_Xref[0] = XX[0]-Xref[0];
    X_Xref[1] = XX[1]-Xref[1];
    X_Xref[2] = XX[2]-Xref[2];
    mul[0] = X_Xref[0]*Q[0]+X_Xref[1]*Q[3]+X_Xref[2]*Q[6];
    mul[1] = X_Xref[0]*Q[1]+X_Xref[1]*Q[4]+X_Xref[2]*Q[7];
    mul[2] = X_Xref[0]*Q[2]+X_Xref[1]*Q[5]+X_Xref[2]*Q[8];
    cost_state = X_Xref[0]*mul[0]+mul[1]*X_Xref[1]+mul[2]*X_Xref[2];
    cost = cost + cost_state;
    B[0] = cos(XX[2]);
    B[1] = -sin(XX[2]);
    B[2] = 0;
    B[3] = sin(XX[2]);
    B[4] = cos(XX[2]);
    B[5] = 0;
    B[6] = 0;
    B[7] = 0;
    B[8] = 1;
    B_U[0] = B[0]*U[3]+B[1]*U[4]+B[2]*U[5];
    B_U[1] = B[3]*U[3]+B[4]*U[4]+B[5]*U[5];
    B_U[2] = B[6]*U[3]+B[7]*U[4]+B[8]*U[5];
    Xdot[0] = B_U[0];
    Xdot[1] = B_U[1];
    Xdot[2] = B_U[2];
    XX[0] = XX[0]+ts*Xdot[0];
    XX[1] = XX[1]+ts*Xdot[1];
    XX[2] = XX[2]+ts*Xdot[2];

    // iteration 03
    X_Xref[0] = XX[0]-Xref[0];
    X_Xref[1] = XX[1]-Xref[1];
    X_Xref[2] = XX[2]-Xref[2];
    mul[0] = X_Xref[0]*Q[0]+X_Xref[1]*Q[3]+X_Xref[2]*Q[6];
    mul[1] = X_Xref[0]*Q[1]+X_Xref[1]*Q[4]+X_Xref[2]*Q[7];
    mul[2] = X_Xref[0]*Q[2]+X_Xref[1]*Q[5]+X_Xref[2]*Q[8];
    cost_state = X_Xref[0]*mul[0]+mul[1]*X_Xref[1]+mul[2]*X_Xref[2];
    cost = cost + cost_state;
    B[0] = cos(XX[2]);
    B[1] = -sin(XX[2]);
    B[2] = 0;
    B[3] = sin(XX[2]);
    B[4] = cos(XX[2]);
    B[5] = 0;
    B[6] = 0;
    B[7] = 0;
    B[8] = 1;
    B_U[0] = B[0]*U[6]+B[1]*U[7]+B[2]*U[8];
    B_U[1] = B[3]*U[6]+B[4]*U[7]+B[5]*U[8];
    B_U[2] = B[6]*U[6]+B[7]*U[7]+B[8]*U[8];
    Xdot[0] = B_U[0];
    Xdot[1] = B_U[1];
    Xdot[2] = B_U[2];
    XX[0] = XX[0]+ts*Xdot[0];
    XX[1] = XX[1]+ts*Xdot[1];
    XX[2] = XX[2]+ts*Xdot[2];

    // iteration 04
    X_Xref[0] = XX[0]-Xref[0];
    X_Xref[1] = XX[1]-Xref[1];
    X_Xref[2] = XX[2]-Xref[2];
    mul[0] = X_Xref[0]*Q[0]+X_Xref[1]*Q[3]+X_Xref[2]*Q[6];
    mul[1] = X_Xref[0]*Q[1]+X_Xref[1]*Q[4]+X_Xref[2]*Q[7];
    mul[2] = X_Xref[0]*Q[2]+X_Xref[1]*Q[5]+X_Xref[2]*Q[8];
    cost_state = X_Xref[0]*mul[0]+mul[1]*X_Xref[1]+mul[2]*X_Xref[2];
    cost = cost + cost_state;
    B[0] = cos(XX[2]);
    B[1] = -sin(XX[2]);
    B[2] = 0;
    B[3] = sin(XX[2]);
    B[4] = cos(XX[2]);
    B[5] = 0;
    B[6] = 0;
    B[7] = 0;
    B[8] = 1;
    B_U[0] = B[0]*U[9]+B[1]*U[10]+B[2]*U[11];
    B_U[1] = B[3]*U[9]+B[4]*U[10]+B[5]*U[11];
    B_U[2] = B[6]*U[9]+B[7]*U[10]+B[8]*U[11];
    Xdot[0] = B_U[0];
    Xdot[1] = B_U[1];
    Xdot[2] = B_U[2];
    XX[0] = XX[0]+ts*Xdot[0];
    XX[1] = XX[1]+ts*Xdot[1];
    XX[2] = XX[2]+ts*Xdot[2];

    // iteration 05
    X_Xref[0] = XX[0]-Xref[0];
    X_Xref[1] = XX[1]-Xref[1];
    X_Xref[2] = XX[2]-Xref[2];
    mul[0] = X_Xref[0]*Q[0]+X_Xref[1]*Q[3]+X_Xref[2]*Q[6];
    mul[1] = X_Xref[0]*Q[1]+X_Xref[1]*Q[4]+X_Xref[2]*Q[7];
    mul[2] = X_Xref[0]*Q[2]+X_Xref[1]*Q[5]+X_Xref[2]*Q[8];
    cost_state = X_Xref[0]*mul[0]+mul[1]*X_Xref[1]+mul[2]*X_Xref[2];
    cost = cost + cost_state;
    B[0] = cos(XX[2]);
    B[1] = -sin(XX[2]);
    B[2] = 0;
    B[3] = sin(XX[2]);
    B[4] = cos(XX[2]);
    B[5] = 0;
    B[6] = 0;
    B[7] = 0;
    B[8] = 1;
    B_U[0] = B[0]*U[12]+B[1]*U[13]+B[2]*U[14];
    B_U[1] = B[3]*U[12]+B[4]*U[13]+B[5]*U[14];
    B_U[2] = B[6]*U[12]+B[7]*U[13]+B[8]*U[14];
    Xdot[0] = B_U[0];
    Xdot[1] = B_U[1];
    Xdot[2] = B_U[2];
    XX[0] = XX[0]+ts*Xdot[0];
    XX[1] = XX[1]+ts*Xdot[1];
    XX[2] = XX[2]+ts*Xdot[2];

  
    return cost;
}






void mpc_controller(double *inputs_vx, double *inputs_vy, double *inputs_omega)
{
  // some initialization iteration states update and nlopt parameters
  unsigned n = 15;                                                                         // number of decision variables
  double lb_var = -0.15;
  double up_var = 0.15;
  double lb[] = {lb_var,lb_var,lb_var,   lb_var,lb_var,lb_var,   lb_var,lb_var,lb_var,   lb_var,lb_var,lb_var,   lb_var,lb_var,lb_var}; // low bound for decision variables
  double ub[] = {up_var,up_var,up_var,   up_var,up_var,up_var,   up_var,up_var,up_var,   up_var,up_var,up_var,   up_var,up_var,up_var};                // up bound for decision variables
  double minf;                                                                             // final optimized J result
  double U[15] = {0,0,0,   0,0,0,   0,0,0,   0,0,0,   0,0,0};
  nlopt_opt opt;

  // main nlopt loop
  opt = nlopt_create(NLOPT_LN_COBYLA, n);                                              // generate the optimization problem using nlopt api
  nlopt_set_lower_bounds(opt, lb);                                                     // set the low bound to optimizer
  nlopt_set_upper_bounds(opt, ub);                                                     // set the up bound to optimizer
  nlopt_set_min_objective(opt, mycost, NULL);                                          // set the cost function
  nlopt_set_xtol_rel(opt, 1e-4);                                                       // set the optimizer valve to converge
  if (nlopt_optimize(opt, U, &minf) < 0) 
  {
      // printf("nlopt failed!\n");
  }
  else 
  {
      printf("we have the optimized min(J) is %f\n", minf);
      printf("in this case we have corresponding decision variables are:\n");
      printf("*****************\n");
      printf("%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n", U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14]);
      printf("*****************\n");
  }
  nlopt_destroy(opt);
  // double vx = U[0];
  // double vy = U[1];
  // double omega = U[2];
  // double mpc_inputs[3] = {vx, vy, omega};
  // return vx;
  *inputs_vx = U[0];
  *inputs_vy = U[1];
  *inputs_omega = U[2];
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //


















































void my_controller(int motiontime,double quat_act[4], double omega_act[3],
                   double q_act[nact], double u_act[nact], double tau[nact],
                  double Kp[nact],double Kd[nact],
                  double q_ref[nact], double u_ref[nact],double tauEst[nact],
                  uint16_t footForce[4], double bodyWorldPos[3]) //, double tracking_error[2])
{
  int i,j;
  double tau_r[3]={0}, kp_r[3]={0}, kd_r[3] = {0};
  double F[4][3]={0};

  double ttime = motiontime*dtime;
  printf("check motiontime %d\n", motiontime);
  printf("check ttime %f\n", ttime);

  //time starts
  //struct timeval begin, end;
  //gettimeofday(&begin, 0);

  //printf("%f %f %f %f \n",quat_act[0],quat_act[1],quat_act[2],quat_act[3]);
  get_stance_forces(quat_act,omega_act,q_act,u_act,F);
  //ram_display(&F[0][0],4,3,"F");

  for (j=0;j<4;j++)
  {
    double qr_ref[3]={0}, ur_ref[3]={0};
    double qr_act[3]={0}, ur_act[3]={0};
    for (i=0;i<3;i++)
    {
      qr_act[i]=q_act[index_leg_act[j][i]];
      ur_act[i]=u_act[index_leg_act[j][i]];
    }

    double force[3] = {0};
    for (i=0;i<3;i++)
      force[i] = F[j][i];


    fsm[j] = leg_state_machine(j,fsm[j],motiontime,quat_act,omega_act,qr_act,ur_act,qr_ref,ur_ref,
                                 tau_r,kp_r,kd_r,footForce[j],force);

    double ts = 4; // > t_free + t_trans + t_stance
    if (ttime >= ts && ttime <= ts+tend) //ensures trajectory is provided after ts seconds
    {
      if (prev_step < step_no) //ensures values are set once per step
      {
        prev_step = step_no;

        double x_ref, y_ref;
        double xdot_ref, ydot_ref;
        double theta_ref, thetadot_ref;


        double x_c, y_c,theta;
        double ang1,ang2;
        double dc[3][3];

        //get trajectory reference
        double x_center = -a;
        double y_center = 0;
        get_trajectory(ttime,ts,tend,&x_ref,&y_ref,
                        &xdot_ref,&ydot_ref,
                        &theta_ref,&thetadot_ref,
                        &theta_ref0,a,x_center,y_center);

        // update the most global reference part for reference information go inside mycost
        x_ref_global = x_ref;
        y_ref_global = y_ref;
        theta0_ref_global = theta_ref;
        // printf("check the reference information %f %f %f\n", x_ref, y_ref, theta_ref);

        // check the reference information
        // printf("ref is %f %f\n", x_ref, y_ref);

        //get orientation (yaw)
        sdquat2dc(quat_act[1],quat_act[2],quat_act[3],quat_act[0],dc); sdprinterr(stderr);
        sddc2ang(dc,&ang1,&ang2,&theta); sdprinterr(stderr);

        if (flag_constant)
          theta0 = theta;
        else
          theta0 += t_stance*omega_act[2];
        // printf("%f %f %f %f\n",theta0,theta_ref,omega_act[2],thetadot_ref);
        double c = cos(theta);
        double s = sin(theta);

        //get world coordinates
        x_c = bodyWorldPos[0];
        y_c = bodyWorldPos[1];

        //get position of point P
        double x_p = x_c + c*px - s*py;
        double y_p = y_c + s*px + c*py;
        // update the most global states variable to let them go inside the mycost function;
        x_p_global = x_p;
        y_p_global = y_p;
        theta0_global = theta0;



























        int flag_mpc = 0;
        if (flag_mpc!=0)
        {
        //Xdot = Xdot_ref + Kp * (X- X_ref) //feedback linearization
        double xdot_p = xdot_ref + Kp_ik*(x_ref - x_p);
        double ydot_p = ydot_ref + Kp_ik*(y_ref - y_p);
        double thetadot_p = thetadot_ref + Kp_ik*(theta_ref - theta0);

        //Xdot = B*U where U = [vx; vy; omega]
        double a11 = c; double a12 = s; double a13 = py;
        double a21 = -s; double a22 =  c; double a23 = -px;
        double a31 = 0;  double a32 = 0; double a33 = 1;

        //U = inv(B)*Xdot
        // Here A = inv(B)

        //MPC: Your vx, vy, omega will come from MPC
        vx =    a11*xdot_p + a12*ydot_p + a13*thetadot_p;
        vy =    a21*xdot_p + a22*ydot_p + a23*thetadot_p;
        omega = a31*xdot_p + a32*ydot_p + a33*thetadot_p;
        // printf("desired inputs from pranav's %f %f %f\n", vx, vy, omega);
       }
       else
       {
         /* basic format from pranav
         //X_P is the measured position
         //Use X_P = [x_p y_p theta] (current state) to solve the MPC here

         //MPC will output vx, vy, omega
         // vx = ? ;
         // vy = ?;
         // omega = ?;
         */

        // - we have measure position [x_p, y_p, theta0]
        double mpc_controls[3];
        double inputs_vx;
        double inputs_vy;
        double inputs_omega;
        mpc_controller(&inputs_vx, &inputs_vy, &inputs_omega);
        vx = inputs_vx;
        vy = inputs_vy;
        omega = inputs_omega;
        // printf("check here %f\n", mpc_controls);


        // vx = 0.05;
        // vy = 0.05;



       }

        //printf("%f %f %f \n",x_ref - x_p, y_ref - y_p, theta_ref-theta);
        fprintf(fidC,"%d %f %f %f\n",step_no, vx, vy, omega);
        fprintf(fidmat,"%d %f %f %f %f %f %f %f %f %f %f %f; \n",
        step_no, x_ref, y_ref, theta_ref,
                           x_p, y_p,
                           x_ref - x_p, y_ref - y_p, theta_ref-theta0,
                           vx, vy, omega);

      }
    }
    
    
    
    
    












































    
    
    else if (ttime > ts+tend)
    {
      if (flag_openmat==1)
      {
        fprintf(fidmat,"];");
        fclose(fidmat);
        flag_openmat = 0;
      }

      if (flag_C==1)
      {
        //fprintf(fidmat,"];");
        fclose(fidC);
        flag_C = 0;
      }

      vx = 0;
      omega = 0;
      vy = 0;
    }



    for (i=0;i<3;i++)
    {
      q_ref[index_leg_act[j][i]] = qr_ref[i];
      u_ref[index_leg_act[j][i]] = ur_ref[i];
      tau[index_leg_act[j][i]] = tau_r[i];
      Kp[index_leg_act[j][i]] = kp_r[i];
      Kd[index_leg_act[j][i]] = kd_r[i];
    }
  }

  //timer ends
  //gettimeofday(&end, 0);
  //long seconds = end.tv_sec - begin.tv_sec;
  //long microseconds = end.tv_usec - begin.tv_usec;
  //double elapsed = seconds + microseconds*1e-6;
  //printf("Time measured: %.3f milli-seconds.\n", elapsed*1e3);


  //printf("%f %f %f \n",q_ref[0],q_ref[1],q_ref[2]);

  //feedback_control(motiontime,q_ref,u_ref,q_act,u_act,tau,Kp,Kd,tauEst);

}


