    // B657 assignment 2 skeleton code
    //
    // Compile with: "make"
    //
    // See assignment handout for command line and project specifications.
    /*
     References:
     [1] http://obsessive-coffee-disorder.com/rgb-to-grayscale-using-cimg/
     [2]Photo new1.jpg by Rachel Moenning on Unsplash
     [3]Photo new2.jpg by Anna Sullivan on Unsplash
    */
    //Link to the header file
    //#define  cimg_use_jpeg
    #include "CImg.h"
    #include <iostream>
    #include <stdlib.h>
    #include <string>
    #include <vector>
    #include <Sift.h>
    #include <sstream>
    #include <math.h>
    #include <cmath>
    #include <time.h>       /* time */
    #include <set>
    //Use the cimg namespace to access functions easily
    using namespace cimg_library;
    using namespace std;

    void draw_descriptor_image(CImg<double> image, const vector<SiftDescriptor> descriptors, const char *filename)
    {
      for(unsigned int i=0; i < descriptors.size(); i++)
        {
          int tx1 = 0, ty1 = 0, tx2 = 0, ty2 = 0;
          double color_point[] = {255.0, 255.0, 0};
          for(int x=-2; x<3; x++)
        for(int y=-2; y<3; y++)
          if(x==0 || y==0)
            for(int c=0; c<3; c++){
              //Find if coordinates are in workspace to draw crosshair
              tx1 = (descriptors[i].col + y - 1);
              ty1 = (descriptors[i].row + x - 1);
              if (tx1 >= 0 && tx1 < image.width() && ty1 >= 0 && ty1 < image.height())
            image( tx1, ty1, 0, c) = color_point[c];                
            }
        }
      image.get_normalize(0,255).save(filename);
    }
    //////////////// Functions for Part 1 start here///////////////////////////////////////////// 
    //Used to calculate new co_ordinate given transformation matrix
    //Formula from Computer Vision: Algorithms and Applications by Richard Szeliski
    void matrix_mul(double const matrix_1[3][3],double const matrix_2[3],double output[2] )
    {
        for(int i=0; i < 2; i ++)
        {
            double denominator=0.0;
            for(int j=0; j<3; j++)
            {
                denominator+=matrix_1[2][j]*matrix_2[j];
                output[i]+=matrix_1[i][j]*matrix_2[j];
            }
            output[i]=output[i]/denominator;
            
        }
    }

    void transpose_array(double a[3][3], double trans[3][3])
    {
        for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j)
            {
                trans[j][i]=a[i][j];
            }
    }

    //source: http://www.cplusplus.com/forum/general/223080/
    void inverse(double const mat[3][3], double inv[3][3])
    {
        //finding determinant
        double determinant=0.0;
        for(int i = 0; i < 3; i++)
            determinant = determinant + (mat[0][i] * (mat[1][(i+1)%3] * mat[2][(i+2)%3] - mat[1][(i+2)%3] * mat[2][(i+1)%3]));
        
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++)
                inv[i][j]=((mat[(j+1)%3][(i+1)%3] * mat[(j+2)%3][(i+2)%3]) - (mat[(j+1)%3][(i+2)%3] * mat[(j+2)%3][(i+1)%3]))/ determinant;
        }
        
    }
    double distance(int x1, int y1,int x2, int y2)
    {
        return (pow(x1-x2,2)+pow(y1-y2,2));
    }

    //
    void find_transformation_matrix(double ip_points[4][2],double result[4][2],double parameters[8])
    {
        CImg<double> X(8,1,1,1,0.0);
        CImg<double> A(8,8,1,1,0.0);
        CImg<double> B(1,8,1,1,0.0);
        
        //Writing the coefficient matrix
        for (int i=0;i<8;i++)
        {
            {
                if(i%2==0)
                { 
                    A(i,0)=ip_points[i/2][0];
                    A(i,1)=ip_points[i/2][1];
                    A(i,2)=1;
                    A(i,6)=-(ip_points[i/2][0]*result[i/2][0]);
                    A(i,7)=-(ip_points[i/2][1]*result[i/2][0]);
                }
                else{
                    A(i,3)=ip_points[i/2][0];
                    A(i,4)=ip_points[i/2][1];
                    A(i,5)=1;
                    A(i,6)=-(ip_points[i/2][0]*result[i/2][1]);
                    A(i,7)=-(ip_points[i/2][1]*result[i/2][1]);
                }
            }
        }
        A.transpose();
        int counter=0;
        for (int i=0;i<4;i++)
        {
            for (int j=0;j<2;j++)
            {
                B(0,counter)=result[i][j];
                counter++;
            }
        }
        
        X=B.solve(A);
        
        for(int i=0;i<8;i++)
        {
            parameters[i]=X(i);
        }
    }

    void image_transformation(CImg<double> input_image,CImg<double> &transformed_img,double matrix[3][3])
    {   

        
        for(int i=0;i<input_image.width();i++)
        {
            for(int j=0; j< input_image.height();j++)
            {
                int x,y;
                double old_pos[]={i,j,1};
                double output[]={0.0,0.0};
                matrix_mul(matrix,old_pos,output);
                x=round(output[0]);
                y=round(output[1]);
                
                if(x>=input_image.width() or y>=input_image.height()or x<=0 or y<=0)
                {}
                else{
                    transformed_img(i,j,0,0)=input_image(x,y,0,0);
                    transformed_img(i,j,0,1)=input_image(x,y,0,1);
                    transformed_img(i,j,0,2)=input_image(x,y,0,2);  
                }
            }
        }

    }
    void image_transformation_billboard(CImg<double> input_image,CImg<double> &transformed_img,double matrix[3][3])
    {
        int count=0;
        
        for(int i=0;i<transformed_img.width();i++)
        {
            for(int j=0; j< transformed_img.height();j++)
            {
                if (transformed_img(i,j,0,0)>150 && transformed_img(i,j,0,1)>150 && transformed_img(i,j,0,2)>150)
                {
                    count++;
                    int x,y;
                    double old_pos[]={i,j,1};
                    double output[]={0.0,0.0};
                    matrix_mul(matrix,old_pos,output);
                    x=round(output[0]);
                    y=round(output[1]);
                    
                    if(x>=input_image.width() or y>=input_image.height()or x<=0 or y<=0)
                    {}
                    else{
                        transformed_img(i,j,0,0)=input_image(x,y,0,0);
                        transformed_img(i,j,0,1)=input_image(x,y,0,1);
                        transformed_img(i,j,0,2)=input_image(x,y,0,2);  
                    }
                }
            }
        }
    }

    //Used this function to find the co-ordinates of the billboard.
    void find_billboard(CImg<double> image, int corners[4][2])
    {
        
        int dist_tl=std::numeric_limits<double>::max();
        int dist_tr=std::numeric_limits<double>::max();
        int dist_bl=std::numeric_limits<double>::max();
        int dist_br=std::numeric_limits<double>::max();
        
        int ltx=0;
        int lty=0;
        int rtx=image.width()-1;
        int rty=0;
        int lbx=0;
        int lby=image.height()-1;
        int rbx=image.width()-1;
        int rby=image.height()-1;
        for(int i=575; i <image.width()-300;i++)
        {
            for(int j=200;j<image.height()-450;j++)
            {
                if (image(i,j,0,0,1)==255 && image(i,j,0,0,2)==255 && image(i,j,0,0,3)==255)
                { 
                    int d1=distance(i,j,ltx,lty);
                    int d2=distance(i,j,rtx,rty);
                    int d3=distance(i,j,lbx,lby);
                    int d4=distance(i,j,rbx,rby);
                    
                    if (d1< dist_tl)
                    {
                        corners[0][0]=i;
                        corners[0][1]=j;
                        dist_tl=d1;
                    }
                    if (d2< dist_tr)
                    {
                        corners[1][0]=i;
                        corners[1][1]=j;
                        dist_tr=d2;
                    }if (d3< dist_bl)
                    {
                        corners[2][0]=i;
                        corners[2][1]=j;
                        dist_bl=d3;
                    }if (d4< dist_br)
                    {
                        corners[3][0]=i;
                        corners[3][1]=j;
                        dist_br=d4;
                    }
                }
                else{
                    image(i,j,0,0)=0;
                    image(i,j,0,1)=0;
                    image(i,j,0,2)=0;
                }
            }
        }
    }

    void to_matrix(double p[8],double h[3][3] )
    {
        int counter=0;
        for (int i=0;i<3;i++)
        {
            for (int j=0;j<3;j++)
            {
                if(i==2 and j==2)
                {
                    h[i][j]=1;  
                }
                else{
                    h[i][j]=p[counter];
                }
                counter++;
            }
        }
        
    }
    /////////////////////////////// Functions for part 1 ends here//////////////////////////////////

    ///////////////////////////////Functions for part 3 starts here/////////////////////////////////

    class sift_match{
    public:
        SiftDescriptor first;
        SiftDescriptor second;
        double dist;
    };
    double sift_distance(SiftDescriptor sd1,SiftDescriptor sd2)
    {
        int sum=0;
        for (int i=0; i<128;i++)
        {
            sum+=pow((sd1.descriptor[i]-sd2.descriptor[i]),2);
        }
        return pow(sum,0.5);
    } 
    vector <sift_match>  get_match(vector<SiftDescriptor> image1_descriptors,vector<SiftDescriptor> image2_descriptors,double threshold){
        vector<sift_match> output;
        
        //First Closest
        SiftDescriptor sd1;
        //Second closest
        SiftDescriptor sd2;
        sift_match temp; 
        
        for (int i=0; i <image1_descriptors.size(); i++)
        {
            double d1=2147483647;
            double d2=2147483647;
            for(int j=0; j<image2_descriptors.size(); j++)
            {
                double distance=sift_distance(image1_descriptors[i],image2_descriptors[j]);
                if(distance<d1)
                {
                    d2=d1;
                    d1=distance;
                    
                    sd2=sd1;
                    sd1=image2_descriptors[j];
                    
                    
                }else if(distance<d2)
                {
                    d2=distance;
                    sd2=image2_descriptors[j];
                }
            }
            
            if(d2==0 || d1/d2<threshold)
            {
                temp.first=image1_descriptors[i];
                temp.second=sd1;
                temp.dist=d1;
                output.push_back(temp);
            }
        }
        return output;
    }

    void pair_images(CImg<double> img1, CImg<double> img2 ,CImg<double> &output_img, vector<sift_match> matches)
    {
        const unsigned char color_line[] = {0, 255, 0};
        
        for (int i =0; i<img1.width();i++)
        {
            for(int j=0 ; j<img1.height();j++)
            {
                output_img(i,j,0,0)=img1(i,j,0,0);
                output_img(i,j,0,1)=img1(i,j,0,1);
                output_img(i,j,0,2)=img1(i,j,0,2);        
            }
        }
        int shift_x=img1.width();
        
        for (int i=0; i < img2.width();i++)
        {
            for (int j=0;j<img2.height();j++)
            {
                output_img(i+shift_x,j,0,0)=img2(i,j,0,0);
                output_img(i+shift_x,j,0,1)=img2(i,j,0,1);
                output_img(i+shift_x,j,0,2)=img2(i,j,0,2);
            }
        }
        for(int i =0;i<matches.size();i++)
        {
            output_img.draw_line(int(matches[i].first.col),int(matches[i].first.row),int(matches[i].second.col)+shift_x,int(matches[i].second.row),color_line);
        }
    }
    // Reference:
    // Random numbers: http://www.cplusplus.com/reference/cstdlib/rand/

    ///////
    void ransac(vector<sift_match> match_points,double best_homography[3][3],vector<sift_match> &ransac_match)
    {
        int MAX_INLIER=-2147483648;
        int counter=0;
        
        //Threshold set to 3 pixels
        int threshold=4; //3 5 9 13 17 34 72;
        for(int iterations=0;iterations<72;iterations++)
        {
            int inlier_count=0;
            std::vector<sift_match>  temp_inliers;
            //generate random 4 matches from the vector of matched points
            double image1[4][2];
            double image2[4][2];
            double transform_matrix[3][3];
            int track_match_index;
            set<int> distinct_indices;
            set<int>::iterator it;
            //srand(int(time(NULL)));
            while(true){
                int random =rand() % match_points.size();
                distinct_indices.insert(random);
                if (distinct_indices.size()==4)
                {
                    break;
                }
            }
            int k=0;
            for (it=distinct_indices.begin(); it!=distinct_indices.end(); ++it)
            {
                track_match_index= *it;
                //cout<<track_match_index<<endl;
                image1[k][0] = match_points[track_match_index].first.col;
                image1[k][1] = match_points[track_match_index].first.row;
                image2[k][0] = match_points[track_match_index].second.col;
                image2[k][1] = match_points[track_match_index].second.row;
                k++;
            }
            double parameters[8];
            find_transformation_matrix(image1,image2,parameters);
            to_matrix(parameters,transform_matrix);
            
            //Test for each remaining sift pair
            //Calculate number of inliers and outliers
            
            for (int i=0; i<match_points.size(); i++)
            {
                
                int img1_x=match_points[i].first.col;
                int img1_y=match_points[i].first.row;
                
                int img2_x=match_points[i].second.col;
                int img2_y=match_points[i].second.row;
                
                double output[2];
                double img1_points[3]={img1_x,img1_y,1};
                
                matrix_mul(transform_matrix,img1_points,output);
                int x_dash=round(output[0]);
                int y_dash=round(output[1]);
                //cout<<x_dash<<" "<<y_dash<<endl;
                
                if (abs(x_dash-img2_x)< threshold && abs(y_dash-img2_y)< threshold)
                {
                    inlier_count++;
                    temp_inliers.push_back(match_points[i]);
                }
                
            }
            //cout<<"countinlier:"<<inlier_count<<" "<<MAX_INLIER<<endl;
            
            if (inlier_count>MAX_INLIER)
            {
                counter++;
                MAX_INLIER=inlier_count;
                ransac_match.clear();
                for (int i = 0; i <temp_inliers.size() ; ++i)
                {
                   
                    ransac_match.push_back(temp_inliers[i]);
                }
                
                for(int i=0;i<3;i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        best_homography[i][j]=transform_matrix[i][j];
                    }
                }
                
            }
            
            
        }
    }
            


    /////////////////////////////// Functions for part 3 ends here//////////////////////////////////


    /////////////////////////////// Functions for part 3 ends here//////////////////////////////////
    /////////Functions for part2 start here////////////
    //downsample by 2
    CImg<double> down_sample_image(CImg<double> image,int channel){
        CImg<double> Down(ceil(image.width()/2.0), ceil(image.height()/2.0), 1, 3, 0);
        for(int i=0;i<ceil(image.width()/2.0);i++){
            for(int j=0;j<ceil(image.height()/2.0);j++){
                
                Down(i,j,0,channel)=image(i*2,j*2,0,channel);
                
            }
        }
        return Down;
    }

    CImg<double> correctDimension(CImg<double> correctimage,int channel){
        CImg<double> correct(correctimage.width()-1, correctimage.height()-1, 1, 3, 0);
        for(int i=0;i<correctimage.width()-1;i++){
            for(int j=0;j<correctimage.height()-1;j++){
                correct(i,j,0,channel)=correctimage(i,j,0,channel);
            }
        }
        return correct;
    }

    //upscale by 2
    CImg<double> Up_scale_image(CImg<double> image,int desired_dim,int channel){
        CImg<double> Up(image.width()*2, image.height()*2, 1, 3, 0);
        for(int i=0;i<image.width();i++){
            for(int j=0;j<image.height();j++){
                Up(i*2,j*2,0,channel)=image(i,j,0,channel);
            }
        }
        if(image.height()*2==desired_dim){
        
            return Up;
        }
        else
        {
            CImg<double> correctedUpscale =correctDimension(Up,channel);
            return correctedUpscale;
        }
        
    }
    /////////////function for part2 end here/////////////////////
    /////part4here////////

    CImg<double> matrixtocimg(double mat[3][3]){
        CImg<double>interimImage(3,3);
    for (int m = 0; m < 3; m++){
        for (int n = 0; n < 3; n++)
        {
            interimImage(m,n)=mat[m][n];
        }
    }

        return interimImage;
        
    }
    void imgtomatrix(CImg<double> image,double mat[3][3]){
        for (int m = 0; m < 3; m++){
            for (int n = 0; n < 3; n++)
            {
                mat[m][n]=image(m,n);
            }
        }
    }

    void image_transformation_p4(CImg<double> input_image,CImg<double> &transformed_img,double matrix[3][3])
    {   

        
        for(int i=0;i<transformed_img.width();i++)
        {
            for(int j=0; j< transformed_img.height();j++)
            {
                int x,y;
                double old_pos[]={i,j,1};
                double output[]={0.0,0.0};
                matrix_mul(matrix,old_pos,output);
                x=round(output[0]);
                y=round(output[1]);
                
                if(x)
                if(x>=input_image.width() or y>=input_image.height()or x<=0 or y<=0)
                {}
                else{
                    transformed_img(i,j,0,0)=input_image(x,y,0,0);
                    transformed_img(i,j,0,1)=input_image(x,y,0,1);
                    transformed_img(i,j,0,2)=input_image(x,y,0,2);
                }
            }
        }

    }

CImg<double> remove_padding(CImg<double> image)
    {
        int top=99999;
        int bottom=00;
        int left=9999999;
        int right=000;


        for(int i=0 ; i<image.width() ; i++)
        {
            for(int j=0 ; j< image.height() ; j++)
            {
                if (image(i,j,0,1)!=0 or image(i,j,0,1)!=0 or image(i,j,0,1)!=0)
                {
                    if(i<left)
                    {
                        left=i;
                    }

                    if(i>right)
                    {
                        right=i;
                    }

                    if(j<top)
                    {
                        top=j;
                    }

                    if(j>bottom)
                    {
                        bottom=j;
                    }
                }      
            }
        }


        CImg<double> temp(right-left,bottom-top,1,3,0);
        for(int i=0 ; i < temp.width() ; i++ )
        {
            for(int j=0; j< temp.height() ; j++)
            {
                temp(i,j,0,0)=image(left+i,top+j,0,0);
                temp(i,j,0,1)=image(left+i,top+j,0,1);
                temp(i,j,0,2)=image(left+i,top+j,0,2);

            }
        }
        return temp;
    }



    ////////////////////////////////Part 4 funtion end here////////////////




    int main(int argc, char **argv)
    {
      try {
        
        /*
          TEST CODE - STARTS
        */
        string part = "";
        /*CImg<double> input_image("images/part2/apple.jpg");
        CImg<double> input_gray = input_image.get_RGBtoHSI().get_channel(2);
        vector<SiftDescriptor> input_descriptors = Sift::compute_sift(input_gray);
        draw_descriptor_image(input_image, input_descriptors, "input_image.jpg");
          */
        /*
          TEST CODE - ENDS
        */
          if(argc < 3)
          {
              cerr << "usage1: " << argv[0] << " part1 " << " poster_input.png " << endl;
              cerr << "usage2: " << argv[0] << " part2 " << " image_1.png "<<" image_2.png "<<" mask.png " << endl;
              cerr << "usage3: " << argv[0] << " part3 " << " image_src.png  "<<" image_dst.png "<< endl;
              cerr << "usage4: " << argv[0] << " part4 " << " image_1.png "<<" image_2.png "<<" image_3.png " << endl;
              return 1;
          }

        part=argv[1];
          
        if(part == "part1"){
          // Billboard
            
            if(argc<3){
                cerr << "usage1: " << argv[0] << " part1 " << " poster_input.png " << endl;
                 return 1;
            }
            string img_path=(argv[2]);
            CImg<double> input_image(("images/part1/" + img_path).c_str());
            CImg<double> input_gray = input_image.get_RGBtoHSI().get_channel(2);
            {
                
                //PART1.1
                CImg<double> lincoln_image("images/part1/lincoln.png");
                CImg<double> transformed_img(lincoln_image.width(),lincoln_image.height(),1,3,0);
                
                //Given Transformation matrix
                double h[3][3]={{0.907,0.258 ,-182},{-0.153,1.44,58.00}, {-0.000306,0.000731,1.0}};
                double h_dash[3][3];
                inverse(h,h_dash);
                
                // Transforming the image given the matrix
                image_transformation(lincoln_image,transformed_img,h_dash);
                transformed_img.save("OutputImages/part1/lincoln_warped.png");
            }
            
            {
                //PART1.2
                double old[4][2]={{318, 256}, {534, 372}, {316, 670}, {73, 473}};
                double new_matrix[4][2]={{141, 131}, {480, 159}, {493, 630},{64, 601}};
                double parameters[8];
                find_transformation_matrix(old,new_matrix,parameters);
                CImg<double> book_image("images/part1/book2.jpg");
                CImg<double> book_result(book_image.width(),book_image.height(),1,3,0);
                
                double h[3][3];
                
                //Initializing matrix for transformation
                cout<<">>>>>>>>>>>>>>Transformation matrix>>>>>>>>>>>>>>>>"<<endl;
                int counter=0;
                for (int i=0;i<3;i++)
                {for (int j=0;j<3;j++)
                {if(i==2 and j==2){
                    h[i][j]=1;
                    cout<<h[i][j]<<"\t";
                }
                else{
                    h[i][j]=parameters[counter];
                    cout<<h[i][j]<<"\t";
                }
                    counter++;
                }
                    cout<<endl;
                }
                image_transformation(book_image,book_result,h);
                book_result.save("OutputImages/part1/book_result.png");
            }
            {
                //PART1.3
                CImg<double> billboard1("images/part1/billboard1.jpg");
                CImg<double> billboard2("images/part1/billboard2.png");
                CImg<double> billboard3("images/part1/billboard3.jpg");
                
                int corners[4][2]={{0,0},{0,0},{0,0},{0,0}};
                double parameters_b1[8];
                double parameters_b2[8];
                double parameters_b3[8];
                
                double old_b1[4][2]={{0,0},{input_image.width()-1,0},{0,input_image.height()-1},{input_image.width()-1,input_image.height()-1}};
                double new_b1[4][2]={{99,60},{535,60},{99,203},{530,202}};
                
                double old_b2[4][2]={{0,0},{input_image.width()-1,0},{0,input_image.height()-1},{input_image.width()-1,input_image.height()-1}};
                double new_b2[4][2]={{177,54},{1105,261},{149,622},{1123,700}};
                
                double old_b3[4][2]={{0,0},{input_image.width()-1,0},{0,input_image.height()-1},{input_image.width()-1,input_image.height()-1}};
                double new_b3[4][2]={{615,285},{1259,260},{610,606},{1259,603}};
                find_transformation_matrix(old_b1,new_b1,parameters_b1);
                find_transformation_matrix(old_b2,new_b2,parameters_b2);
                find_transformation_matrix(old_b3,new_b3,parameters_b3);
                
                double h1[3][3];
                double h2[3][3];
                double h3[3][3];
                
                to_matrix(parameters_b1,h1);
                to_matrix(parameters_b2,h2);
                to_matrix(parameters_b3,h3);
                
                //Initializing matrix for transformation
                double h1_dash[3][3];
                double h2_dash[3][3];
                double h3_dash[3][3];
                
                inverse(h1,h1_dash);
                inverse(h2,h2_dash);
                inverse(h3,h3_dash);
                
                image_transformation_billboard(input_image,billboard1,h1_dash);
                image_transformation_billboard(input_image,billboard2,h2_dash);
                image_transformation_billboard(input_image,billboard3,h3_dash);
                
                billboard1.save("OutputImages/part1/synthetic_billboard1.png");
                billboard2.save("OutputImages/part1/synthetic_billboard2.png");
                billboard3.save("OutputImages/part1/synthetic_billboard3.png");
            }
        }
        else if(part == "part2"){
            if(argc != 5)
            {
                
                cerr << "usage2: " << argv[0] << " part2 " << " image_1.png "<<" image_2.png "<<" mask.png " << endl;
                return 1;
            }

          // Blending
          string Image_1=argv[2];
          string Image_2=argv[3];
          string mask=argv[4];
            
            //initializing RGB vectors for Image1
            CImg<double> input_image1(("images/part2/"+Image_1).c_str()),r1(input_image1.width(), input_image1.height(), 1, 3, 0), g1(input_image1.width(), input_image1.height(), 1, 3, 0), b1(input_image1.width(), input_image1.height(), 1, 3, 0),original1(input_image1.width(), input_image1.height(), 1, 3, 0);
            
            //initializing RGB vectors for Image2
            CImg<double> input_image2(("images/part2/"+Image_2).c_str()),r2(input_image2.width(), input_image2.height(), 1, 3, 0), g2(input_image2.width(), input_image2.height(), 1, 3, 0), b2(input_image2.width(), input_image2.height(), 1, 3, 0),original2(input_image2.width(), input_image2.height(), 1, 3, 0);
           
            //initializing RGB vectors for Mask
            CImg<double> imagemask(("images/part2/"+mask).c_str()),rm(imagemask.width(), imagemask.height(), 1, 3, 0), gm(imagemask.width(), imagemask.height(), 1, 3, 0), bm(imagemask.width(), imagemask.height(), 1, 3, 0);
            
            
            //splitting the image1  into RGB channels
            cimg_forXY(input_image1,x,y){
            r1(x,y,0,0) = input_image1(x,y,0,0),    // Red component of input_image1 sent to r1
            g1(x,y,0,1) = input_image1(x,y,0,1),    // Green component of input_image1 sent to g1
            b1(x,y,0,2) = input_image1(x,y,0,2);    // Blue component of input_image1 sent to b1
            }
            
            //splitting the image2  into RGB channels
            cimg_forXY(input_image2,x,y){
                r2(x,y,0,0) = input_image2(x,y,0,0),    // Red component of input_image2 sent to r1
                g2(x,y,0,1) = input_image2(x,y,0,1),    // Green component of input_image2 sent to g1
                b2(x,y,0,2) = input_image2(x,y,0,2);    // Blue component of input_image2 sent to b1
            }
            
            //splitting the mask  into RGB channels
            for (int x=0;x<imagemask.width();x++){
                for (int y=0;y<imagemask.height();y++){
                        rm(x,y,0,0) = imagemask(x,y),    // component of mask sent to rm
                        gm(x,y,0,1) = imagemask(x,y),    // component of mask sent to gm
                        bm(x,y,0,2) = imagemask(x,y);    // component of mask sent to bm
                }
            }
           //MaskCreator
            /*CImg<double> input_image1_1("new1.jpg");
             CImg<double> input_image1_2("new2.jpg");
             CImg<double> originalm(input_image1_1.width(), input_image1_1.height(), 1, 3, 0);
             
            for (int x=0;x<input_image1_1.width();x++){
                for (int y=0;y<input_image1_1.height();y++){
                    if(y<=input_image1_1.height()/2){
                    originalm(x,y,0,0) = 255,    // component of mask sent to rm
                    originalm(x,y,0,1) = 255,    // component of mask sent to gm
                    originalm(x,y,0,2) = 255;// component of mask sent to bm
                    }
                    else{
                        originalm(x,y,0,0) = 0,    // component of mask sent to rm
                        originalm(x,y,0,1) = 0,    // component of mask sent to gm
                        originalm(x,y,0,2) = 0;// component of mask sent to bm
                    }
                }
            }
            originalm.save("mask2.jpg");*/
            //define gaussian kernel
            CImg<double> GuassianKernel(5,5,1,1,0);
            GuassianKernel(0,0,0,0)=1/256.0;
            GuassianKernel(0,1,0,0)=4/256.0;
            GuassianKernel(0,2,0,0)=6/256.0;
            GuassianKernel(0,3,0,0)=4/256.0;
            GuassianKernel(0,4,0,0)=1/256.0;
            GuassianKernel(1,0,0,0)=4/256.0;
            GuassianKernel(1,1,0,0)=16/256.0;
            GuassianKernel(1,2,0,0)=24/256.0;
            GuassianKernel(1,3,0,0)=16/256.0;
            GuassianKernel(1,4,0,0)=4/256.0;
            GuassianKernel(2,0,0,0)=6/256.0;
            GuassianKernel(2,1,0,0)=24/256.0;
            GuassianKernel(2,2,0,0)=36/256.0;
            GuassianKernel(2,3,0,0)=24/256.0;
            GuassianKernel(2,4,0,0)=6/256.0;
            GuassianKernel(3,0,0,0)=4/256.0;
            GuassianKernel(3,1,0,0)=16/256.0;
            GuassianKernel(3,2,0,0)=24/256.0;
            GuassianKernel(3,3,0,0)=16/256.0;
            GuassianKernel(3,4,0,0)=4/256.0;
            GuassianKernel(4,0,0,0)=1/256.0;
            GuassianKernel(4,1,0,0)=4/256.0;
            GuassianKernel(4,2,0,0)=6/256.0;
            GuassianKernel(4,3,0,0)=4/256.0;
            GuassianKernel(4,4,0,0)=1/256.0;
            
            //Creating Guassian Pyramid for Image1
            int depth=6;
            CImgList< double > GuassianRed1(depth);
            CImgList< double > GuassianGreen1(depth);
            CImgList< double > GuassianBlue1(depth);
            GuassianRed1.insert(r1,0);
            GuassianGreen1.insert(g1,0);
            GuassianBlue1.insert(b1,0);
            
            CImg<double> GuassPy1(GuassianRed1.at(0).width(), GuassianRed1.at(0).height(), 1, 3, 0);
            cimg_forXY(GuassianRed1.at(0),x,y){
                
                GuassPy1(x,y,0,0) = GuassianRed1.at(0)(x,y,0,0);
                GuassPy1(x,y,0,1) = GuassianGreen1.at(0)(x,y,0,1);
                GuassPy1(x,y,0,2) = GuassianBlue1.at(0)(x,y,0,2);
            }
            
            GuassPy1.get_normalize(0,255).save("OutputImages/part2/GuassianPyramidImage1level0.jpg");
            
            for(int i=1;i<depth;i++){
                std::ostringstream s;
                s << i;
                std::string i_as_string=(s.str());
                r1.convolve(GuassianKernel);
                CImg<double> Down1=down_sample_image(r1,0);
                GuassianRed1.insert(Down1,i);
                //Down1.save(("images/part2/G"+i_as_string+"forRed1.jpg").c_str());
                r1=Down1;
                g1.convolve(GuassianKernel);
                CImg<double> Down2=down_sample_image(g1,1);
                GuassianGreen1.insert(Down2,i);
                //Down2.save(("images/part2/G"+i_as_string+"forGreen1.jpg").c_str());
                g1=Down2;
                
                b1.convolve(GuassianKernel);
                CImg<double> Down3=down_sample_image(b1,2);
                GuassianBlue1.insert(Down3,i);
                //Down3.save(("images/part2/G"+i_as_string+"forBlue1.jpg").c_str());
                b1=Down3;
                CImg<double> GuassPy1(GuassianRed1.at(i).width(), GuassianRed1.at(i).height(), 1, 3, 0);
                cimg_forXY(GuassianRed1.at(i),x,y){
                    
                    GuassPy1(x,y,0,0) = GuassianRed1.at(i)(x,y,0,0);
                    GuassPy1(x,y,0,1) = GuassianGreen1.at(i)(x,y,0,1);
                    GuassPy1(x,y,0,2) = GuassianBlue1.at(i)(x,y,0,2);
                }
                
                GuassPy1.get_normalize(0,255).save(("OutputImages/part2/GuassianPyramidImage1level"+i_as_string+".jpg").c_str());
            }
            
            rm/=255.0;
            gm/=255.0;
            bm/=255.0;
            
            //Creating Guassian Pyramid for mask
            CImgList< double > GuassianRedmask(depth);
            CImgList< double > GuassianGreenmask(depth);
            CImgList< double > GuassianBluemask(depth);
            GuassianRedmask.insert(rm,0);
            GuassianGreenmask.insert(gm,0);
            GuassianBluemask.insert(bm,0);
            
            CImg<double> GuassPym(GuassianRedmask.at(0).width(), GuassianRedmask.at(0).height(), 1, 3, 0);
            cimg_forXY(GuassianRedmask.at(0),x,y){
                
                GuassPym(x,y,0,0) = GuassianRedmask.at(0)(x,y,0,0);
                GuassPym(x,y,0,1) = GuassianGreenmask.at(0)(x,y,0,1);
                GuassPym(x,y,0,2) = GuassianBluemask.at(0)(x,y,0,2);
            }
            GuassPym.get_normalize(0,255).save("OutputImages/part2/GuassianPyramidMasklevel0.jpg");
            
            for(int i=1;i<depth;i++){
                std::ostringstream s;
                s << i;
                std::string i_as_string=(s.str());
                rm.convolve(GuassianKernel);
                CImg<double> Down1=down_sample_image(rm,0);
                GuassianRedmask.insert(Down1,i);
                //rm.save(("images/part2/G"+i_as_string+"forRedmask.jpg").c_str());
                rm=Down1;
                gm.convolve(GuassianKernel);
                CImg<double> Down2=down_sample_image(gm,1);
                GuassianGreenmask.insert(Down2,i);
                //Down2.save(("images/part2/G"+i_as_string+"forGreenmask.jpg").c_str());
                gm=Down2;
                bm.convolve(GuassianKernel);
                CImg<double> Down3=down_sample_image(bm,2);
                GuassianBluemask.insert(Down3,i);
                //Down3.save(("images/part2/G"+i_as_string+"forBluemask.jpg").c_str());
                bm=Down3;
                CImg<double> GuassPym(GuassianBluemask.at(i).width(), GuassianBluemask.at(i).height(), 1, 3, 0);
                cimg_forXY(GuassianBluemask.at(i),x,y){
                    
                    GuassPym(x,y,0,0) = GuassianRedmask.at(i)(x,y,0,0);
                    GuassPym(x,y,0,1) = GuassianGreenmask.at(i)(x,y,0,1);
                    GuassPym(x,y,0,2) = GuassianBluemask.at(i)(x,y,0,2);
                }
                
                GuassPym.get_normalize(0,255).save(("OutputImages/part2/GuassianPyramidMasklevel"+i_as_string+".jpg").c_str());
            }
            
            
            //Creating Guassian Pyramid for Image2
            CImgList< double > GuassianRed2(depth);
            CImgList< double > GuassianGreen2(depth);
            CImgList< double > GuassianBlue2(depth);
            GuassianRed2.insert(r2,0);
            GuassianGreen2.insert(g2,0);
            GuassianBlue2.insert(b2,0);
            
            CImg<double> GuassPy2(GuassianRed1.at(0).width(), GuassianRed1.at(0).height(), 1, 3, 0);
            cimg_forXY(GuassianRed2.at(0),x,y){
                
                GuassPy2(x,y,0,0) = GuassianRed2.at(0)(x,y,0,0);
                GuassPy2(x,y,0,1) = GuassianGreen2.at(0)(x,y,0,1);
                GuassPy2(x,y,0,2) = GuassianBlue2.at(0)(x,y,0,2);
            }
            GuassPy2.get_normalize(0,255).save("OutputImages/part2/GuassianPyramidImage2level0.jpg");
            
            for(int i=1;i<depth;i++){
                std::ostringstream s;
                s << i;
                std::string i_as_string=(s.str());
                r2.convolve(GuassianKernel);
                CImg<double> Down1=down_sample_image(r2,0);
                GuassianRed2.insert(Down1,i);
                //r2.save(("images/part2/G"+i_as_string+"forRed2.jpg").c_str());
                r2=Down1;
                g2.convolve(GuassianKernel);
                CImg<double> Down2=down_sample_image(g2,1);
                GuassianGreen2.insert(Down2,i);
                //Down2.save(("images/part2/G"+i_as_string+"forGreen2.jpg").c_str());
                g2=Down2;
                b2.convolve(GuassianKernel);
                CImg<double> Down3=down_sample_image(b2,2);
                GuassianBlue2.insert(Down3,i);
                //Down3.save(("images/part2/G"+i_as_string+"forBlue2.jpg").c_str());
                b2=Down3;
                CImg<double> GuassPy2(GuassianRed2.at(i).width(), GuassianRed2.at(i).height(), 1, 3, 0);
                cimg_forXY(GuassianRed2.at(i),x,y){
                    
                    GuassPy2(x,y,0,0) = GuassianRed2.at(i)(x,y,0,0);
                    GuassPy2(x,y,0,1) = GuassianGreen2.at(i)(x,y,0,1);
                    GuassPy2(x,y,0,2) = GuassianBlue2.at(i)(x,y,0,2);
                }
                
                GuassPy2.get_normalize(0,255).save(("OutputImages/part2/GuassianPyramidImage2level"+i_as_string+".jpg").c_str());
            }
            
            //Create Laplacian pyramid for image1
            
            CImgList< double > LaplacianRed1(depth);
            CImgList< double > LaplacianGreen1(depth);
            CImgList< double > LaplacianBlue1(depth);
            LaplacianRed1.insert(GuassianRed1.at(depth-1),0);
            LaplacianGreen1.insert(GuassianGreen1.at(depth-1),0);
            LaplacianBlue1.insert(GuassianBlue1.at(depth-1),0);
            CImg<double> LapPy1(LaplacianRed1.at(0).width(), LaplacianRed1.at(0).height(), 1, 3, 0);
            cimg_forXY(LaplacianRed1.at(0),x,y){
                
                LapPy1(x,y,0,0) = LaplacianRed1.at(0)(x,y,0,0);
                LapPy1(x,y,0,1) = LaplacianGreen1.at(0)(x,y,0,1);
                LapPy1(x,y,0,2) = LaplacianBlue1.at(0)(x,y,0,2);
            }
            
            LapPy1.save("OutputImages/part2/LaplacianPyramidImage1level0.jpg");
           
            for(int i=1;i<depth;i++){
                std::ostringstream s;
                s << i;
                std::string i_as_string=(s.str());
                CImg<double> upperLevelR1=GuassianRed1.at(depth-i-1);
                CImg<double> lowerLevelR1=4*Up_scale_image(GuassianRed1.at(depth-i),upperLevelR1.height(),0).convolve(GuassianKernel);
                CImg<double> LR1(upperLevelR1.width(),upperLevelR1.height(),1,3,0);
                cimg_forXY(upperLevelR1,x,y){
                    
                    LR1(x,y,0,0) = upperLevelR1(x,y,0,0)-lowerLevelR1(x,y,0,0);
                    
                }
                
                LaplacianRed1.insert(LR1,i);
                //LR1.save(("images/part2/L"+i_as_string+"forRed1.jpg").c_str());
                CImg<double> upperLevelG1=GuassianGreen1.at(depth-i-1);
                CImg<double> lowerLevelG1=4*Up_scale_image(GuassianGreen1.at(depth-i),upperLevelG1.height(),1).convolve(GuassianKernel);
                CImg<double> LG1(upperLevelG1.width(),upperLevelG1.height(),1,3,0);
                cimg_forXY(upperLevelG1,x,y){
                    
                    LG1(x,y,0,1) = upperLevelG1(x,y,0,1)-lowerLevelG1(x,y,0,1);
                    
                }
                LaplacianGreen1.insert(LG1,i);
                //LG1.save(("images/part2/L"+i_as_string+"forGreen1.jpg").c_str());
                CImg<double> upperLevelB1=GuassianBlue1.at(depth-i-1);
                CImg<double> lowerLevelB1=4*Up_scale_image(GuassianBlue1.at(depth-i),upperLevelB1.height(),2).convolve(GuassianKernel);
                CImg<double> LB1(upperLevelB1.width(),upperLevelB1.height(),1,3,0);
                cimg_forXY(upperLevelB1,x,y){
                    
                    LB1(x,y,0,2) = upperLevelB1(x,y,0,2)-lowerLevelB1(x,y,0,2);
                    
                }
                LaplacianBlue1.insert(LB1,i);
                //LB1.save(("images/part2/L"+i_as_string+"forBlue1.jpg").c_str());
                CImg<double> LapPy1(LaplacianRed1.at(i).width(), LaplacianRed1.at(i).height(), 1, 3, 0);
                cimg_forXY(LaplacianRed1.at(i),x,y){
                    
                    LapPy1(x,y,0,0) = LaplacianRed1.at(i)(x,y,0,0);
                    LapPy1(x,y,0,1) = LaplacianGreen1.at(i)(x,y,0,1);
                    LapPy1(x,y,0,2) = LaplacianBlue1.at(i)(x,y,0,2);
                }
                
                LapPy1.get_normalize(0,255).save(("OutputImages/part2/LaplacianPyramidImage1level"+i_as_string+".jpg").c_str());
            }
            
            //Create Laplacian pyramid for image2
            CImgList< double > LaplacianRed2(depth);
            CImgList< double > LaplacianGreen2(depth);
            CImgList< double > LaplacianBlue2(depth);
            LaplacianRed2.insert(GuassianRed2.at(depth-1),0);
            LaplacianGreen2.insert(GuassianGreen2.at(depth-1),0);
            LaplacianBlue2.insert(GuassianBlue2.at(depth-1),0);
            CImg<double> LapPy2(LaplacianRed2.at(0).width(), LaplacianRed2.at(0).height(), 1, 3, 0);
            cimg_forXY(LaplacianRed2.at(0),x,y){
                
                LapPy2(x,y,0,0) = LaplacianRed2.at(0)(x,y,0,0);
                LapPy2(x,y,0,1) = LaplacianGreen2.at(0)(x,y,0,1);
                LapPy2(x,y,0,2) = LaplacianBlue2.at(0)(x,y,0,2);
            }
            
            LapPy2.get_normalize(0,255).save("OutputImages/part2/LaplacianPyramidImage2level0.jpg");
            
             for(int i=1;i<depth;i++){
             std::ostringstream s;
             s << i;
             std::string i_as_string=(s.str());
             CImg<double> upperLevelR2=GuassianRed2.at(depth-i-1);
             CImg<double> lowerLevelR2=4*Up_scale_image(GuassianRed2.at(depth-i),upperLevelR2.height(),0).convolve(GuassianKernel);
             CImg<double> LR2(upperLevelR2.width(),upperLevelR2.height(),1,3,0);
             cimg_forXY(upperLevelR2,x,y){
             
             LR2(x,y,0,0) = upperLevelR2(x,y,0,0)-lowerLevelR2(x,y,0,0);
             
             }
             
             LaplacianRed2.insert(LR2,i);
             //LR1.save(("images/part2/L"+i_as_string+"forRed1.jpg").c_str());
             CImg<double> upperLevelG2=GuassianGreen2.at(depth-i-1);
             CImg<double> lowerLevelG2=4*Up_scale_image(GuassianGreen2.at(depth-i),upperLevelG2.height(),1).convolve(GuassianKernel);
             CImg<double> LG2(upperLevelG2.width(),upperLevelG2.height(),1,3,0);
             cimg_forXY(upperLevelG2,x,y){
             
             LG2(x,y,0,1) = upperLevelG2(x,y,0,1)-lowerLevelG2(x,y,0,1);
             
             }
             LaplacianGreen2.insert(LG2,i);
             //LG1.save(("images/part2/L"+i_as_string+"forGreen1.jpg").c_str());
             CImg<double> upperLevelB2=GuassianBlue2.at(depth-i-1);
             CImg<double> lowerLevelB2=4*Up_scale_image(GuassianBlue2.at(depth-i),upperLevelB2.height(),2).convolve(GuassianKernel);
             CImg<double> LB2(upperLevelB2.width(),upperLevelB2.height(),1,3,0);
             cimg_forXY(upperLevelB2,x,y){
             
             LB2(x,y,0,2) = upperLevelB2(x,y,0,2)-lowerLevelB2(x,y,0,2);
             
             }
             LaplacianBlue2.insert(LB2,i);
             //LB1.save(("images/part2/L"+i_as_string+"forBlue1.jpg").c_str());
             CImg<double> LapPy2(LaplacianRed2.at(i).width(), LaplacianRed2.at(i).height(), 1, 3, 0);
             cimg_forXY(LaplacianRed2.at(i),x,y){
             
             LapPy2(x,y,0,0) = LaplacianRed2.at(i)(x,y,0,0);
             LapPy2(x,y,0,1) = LaplacianGreen2.at(i)(x,y,0,1);
             LapPy2(x,y,0,2) = LaplacianBlue2.at(i)(x,y,0,2);
             }
             
             LapPy2.get_normalize(0,255).save(("OutputImages/part2/LaplacianPyramidImage2level"+i_as_string+".jpg").c_str());
             }
            
            //create laplacian pyramid for blended image
            CImgList< double > LaplacianblendR(depth);
            CImgList< double > LaplacianblendG(depth);
            CImgList< double > LaplacianblendB(depth);
           for(int i=0;i<depth;i++){
               
               CImg<double> Lap_R1=LaplacianRed1.at(i);
               CImg<double> Lap_R2=LaplacianRed2.at(i);
               CImg<double> GM_R=GuassianRedmask.at(depth-i-1);
               CImg<double> Lap_blendR(GM_R.width(),GM_R.height(),1, 3, 0);
               CImg<double> Lap_G1=LaplacianGreen1.at(i);
               CImg<double> Lap_G2=LaplacianGreen2.at(i);
               CImg<double> GM_G=GuassianGreenmask.at(depth-i-1);
               CImg<double> Lap_blendG(GM_G.width(),GM_G.height(),1, 3, 0);
               CImg<double> Lap_B1=LaplacianBlue1.at(i);
               CImg<double> Lap_B2=LaplacianBlue2.at(i);
               CImg<double> GM_B=GuassianBluemask.at(depth-i-1);
               CImg<double> Lap_blendB(GM_B.width(),GM_B.height(),1, 3, 0);
               cimg_forXY(Lap_R1,x,y){
                   Lap_blendR(x,y,0,0)=Lap_R2(x,y,0,0)*GM_R(x,y,0,0)+(1-GM_R(x,y,0,0))*Lap_R1(x,y,0,0);
                   Lap_blendG(x,y,0,1)=Lap_G2(x,y,0,1)*GM_G(x,y,0,1)+(1-GM_G(x,y,0,1))*Lap_G1(x,y,0,1);
                   Lap_blendB(x,y,0,2)=Lap_B2(x,y,0,2)*GM_B(x,y,0,2)+(1-GM_B(x,y,0,2))*Lap_B1(x,y,0,2);
               }
               LaplacianblendR.insert(Lap_blendR,i);
               LaplacianblendG.insert(Lap_blendG,i);
               LaplacianblendB.insert(Lap_blendB,i);
               
           }
            for(int i=0; i< depth ;i++){
                CImg<double> LapBlend(LaplacianblendR.at(i).width(), LaplacianblendR.at(i).height(), 1, 3, 0);
                std::ostringstream s;
                s << i;
                std::string i_as_string=(s.str());
            cimg_forXY(LaplacianblendR.at(i),x,y){
                
                LapBlend(x,y,0,0) = LaplacianblendR.at(i)(x,y,0,0);
                LapBlend(x,y,0,1) = LaplacianblendG.at(i)(x,y,0,1);
                LapBlend(x,y,0,2) = LaplacianblendB.at(i)(x,y,0,2);
            }
            
            LapBlend.get_normalize(0,255).save(("OutputImages/part2/lapblend"+i_as_string+".jpg").c_str());
            }
            
            //Final Blended Image Reconstruction
            for(int i=0; i< depth-1 ;i++){
                LaplacianblendR.at(i+1)=4*Up_scale_image(LaplacianblendR.at(i),LaplacianblendR.at(i+1).height(),0).convolve(GuassianKernel)+LaplacianblendR.at(i+1);
                LaplacianblendG.at(i+1)=4*Up_scale_image(LaplacianblendG.at(i),LaplacianblendG.at(i+1).height(),1).convolve(GuassianKernel)+LaplacianblendG.at(i+1);
                LaplacianblendB.at(i+1)=4*Up_scale_image(LaplacianblendB.at(i),LaplacianblendB.at(i+1).height(),2).convolve(GuassianKernel)+LaplacianblendB.at(i+1);
            }
            
            cimg_forXY(LaplacianblendB.at(depth-1),x,y){
                
                original1(x,y,0,0) = LaplacianblendR.at(depth-1)(x,y,0,0);
                original1(x,y,0,1) = LaplacianblendG.at(depth-1)(x,y,0,1);
                original1(x,y,0,2) = LaplacianblendB.at(depth-1)(x,y,0,2);
            }
            
            original1.get_normalize(0,255).save("OutputImages/part2/FinalBlend.jpg");
            
            
        }
        else if(part == "part3"){
            if(argc != 4){
                
                cerr << "usage3: " << argv[0] << " part3 " << " image_src.png  "<<" image_dst.png "<< endl;
                
                return 1;
            }
          // RANSAC
            string Image1=argv[2];
            string Image2=argv[3];
            //string img_path=(("images/part3/"+Image_1).c_str());
            CImg<double> image_1(("images/part3/" + Image1).c_str());
            CImg<double> image1_gray = image_1.get_RGBtoHSI().get_channel(2);
            
            vector<SiftDescriptor> image1_descriptors = Sift::compute_sift(image1_gray);
            
            //string image_two=(("images/part3/"+Image_2).c_str());
            CImg <double> image_2(("images/part3/" + Image2).c_str());
            CImg<double> image2_gray=image_2.get_RGBtoHSI().get_channel(2);
            
            vector<SiftDescriptor> image2_descriptors = Sift::compute_sift(image2_gray);
            
            //Getting matches for given points
            vector<sift_match> matches=get_match(image1_descriptors,image2_descriptors,0.7);
            
            CImg<double> sift_img(image_1.width()+image_2.width(),max(image_1.height(),image_2.height()),1,3);
            
            pair_images(image_1,image_2,sift_img,matches);
            sift_img.save("OutputImages/part3/sift_matches.jpg");
            //RANSAC PART
            CImg<double> ransac_img(image_1.width()+image_2.width(),max(image_1.height(),image_2.height()),1,3);
            vector<sift_match> ransac_match;
            double h[3][3];
            ransac(matches,h,ransac_match);
            pair_images(image_1,image_2,ransac_img,ransac_match);
            ransac_img.save("OutputImages/part3/ransac_matches.png");
            
            
        }

        else if(part == "part4")
        {
            // Panorama
            
            if(argc<5){
                cerr << "usage4: " << argv[0] << " part4 " << " image_1.png "<<" image_2.png "<<" image_3.png " << endl;
                return 1;}
            std::vector<CImg<double> > images;
            std::vector<vector<SiftDescriptor> >  sift_points;
            std::vector<vector<sift_match> >  sift_matches;
            CImgList<double>Homolist(sift_matches.size());
            CImg<double> latest_img;

            //std::vector<double *[3] > homographies;
            for (int i = 2; i < argc; ++i)
            {
                int center=(argc-((argc-2)/2)-1);
                string img_path=(argv[i]);
                CImg<double> image_1(("images/part4/" + img_path).c_str());
                
                if(i==center)
                {
                    CImg<double> paddedImage(2000,2000,1,3,0);
                    cimg_forXY(paddedImage,x,y)
                    {
                        paddedImage(x,y)=0;
                    }
                    for(int j=0;j<(image_1.width());j++)
                    {
                        for(int k=0; k< image_1.height();k++)
                        {
                            paddedImage(1000-(image_1.width()/2)+j,1000-(image_1.height()/2)+k,0,0)=image_1(j,k,0,0);
                            paddedImage(1000-(image_1.width()/2)+j,1000-(image_1.height()/2)+k,0,1)=image_1(j,k,0,1);
                            paddedImage(1000-(image_1.width()/2)+j,1000-(image_1.height()/2)+k,0,2)=image_1(j,k,0,2);
                        }
                    }

                    latest_img=paddedImage;
                    //latest_img.save("latest_img_test.png");
                    images.push_back(paddedImage);
                }
                else
                {
                    images.push_back(image_1);
                }
                CImg<double> image1_gray = images[i-2].get_RGBtoHSI().get_channel(2);
                vector<SiftDescriptor> image1_descriptors = Sift::compute_sift(image1_gray);
                sift_points.push_back(image1_descriptors);
                
            }
            //images[0].save("Part4Trial.png");
            
            
            int center=sift_points.size()/2;
         
            // From Center to Right
            for (int i=center+1; i <sift_points.size();i++)
            {

                CImg<double> image1_gray = latest_img.get_RGBtoHSI().get_channel(2);
                vector<SiftDescriptor> image1_descriptors = Sift::compute_sift(image1_gray);

                CImg<double> image2_gray = images[i].get_RGBtoHSI().get_channel(2);
                vector<SiftDescriptor> image2_descriptors = Sift::compute_sift(image2_gray);

                double a[3][3];
                vector<sift_match> ransac_match;

                std::vector<sift_match> s=get_match(image2_descriptors,image1_descriptors,0.7);
                ransac(s,a,ransac_match);
                //imgtomatrix(Homolist.at(0),a);
                double a_dash[3][3];
                inverse(a,a_dash);
                //CImg<double> transformed_img(latest_img.width(),latest_img.height(),1,3,0);

                image_transformation_p4(images[i],latest_img,a_dash);
            }

            // From Center to Left
            for (int i=center-1; i >=0;i--)
            {

                CImg<double> image1_gray = latest_img.get_RGBtoHSI().get_channel(2);
                vector<SiftDescriptor> image1_descriptors = Sift::compute_sift(image1_gray);

                CImg<double> image2_gray = images[i].get_RGBtoHSI().get_channel(2);
                vector<SiftDescriptor> image2_descriptors = Sift::compute_sift(image2_gray);

                double a[3][3];
                vector<sift_match> ransac_match;
                std::vector<sift_match> s=get_match(image2_descriptors,image1_descriptors,0.7);
                ransac(s,a,ransac_match);
                //imgtomatrix(Homolist.at(0),a);
                double a_dash[3][3];
                inverse(a,a_dash);
                //CImg<double> transformed_img(latest_img.width(),latest_img.height(),1,3,0);

                image_transformation_p4(images[i],latest_img,a_dash);
                

            }
            CImg<double> clipped_img=remove_padding(latest_img);
            clipped_img.save("OutputImages/part4/panorama.png");
        }   
            
        // feel free to add more conditions for other parts (e.g. more specific)
        //  parts, for debugging, etc.
      }
      catch(const string &err) {
        cerr << "Error: " << err << endl;
      }
    }
