/* Helix Extrude for OpenSCAD
 * Abe Kabakoff (occamsshavingkit@gmail.com)
 * Released 2016-12-29
 * The helix_extrude function takes a shape (a set of 3D points) and twists it up a helix.
 * The user is responsible for calculating the points for the shape he/she wants to use. 
 * he_translate, he_rotate and he_scale are available for doing those functions on a set of 3D points.
 * a multimatrix transform can be done by calling transform_shape
 * A check is performed to ensure the points are on the same plane, but it only issues a warning.
 * There is also he_circle and he_square for the lazy among you. (he_circle obeys $fn, $fa and $fs rules)
 * oh yeah... and
 ***************** THREADS!! ****************
 * use iso_bolt().
 * give it a diameter, a length, a pitch and whether or not it is an internal thread.
 * it will return a ISO profile bolt.
 * if you want a different bolt, make a profile and use helix_extrude.
 * 
 * Functions:
 * helix_extrude(shape, pitch, rotations, convexity) - extrudes a shape
 * transform_shape(matrix, shape) - returns a shape transformed by the matrix
 * he_rotate(rv, shape) - returns a rotated shape. Rotations are applied in the order X, Y, Z.
 * he_translate(tv, shape) - returns a translated shape.
 * he_scale(sv, shape) - returns a scaled shape.
 * he_circle(diameter) - a circle centered on the origin.
 * he_square(v, center) - a square with the size in the vector, centered on the origin or not.
 * matrix_power(matrix, power) - returns a matrix taken to a power.
 * iso_bolt(diameter, pitch, length, internal) - returns a bolt with threads of nominal diameter, with a given pitch, that is a given length.
 *      If internal is true it makes the bolt a bit bigger and it is for subtracting from your solid (using difference()).
 * iso_thread_profile - a helper function for iso_bolt. Not meant for external use.
 */

/* 
 height amount to change on internal threads to give some wiggle room.
 */
internal_fudge = 1/8;

/* takes a matrix to a power
   works only for positive integers!
*/
function matrix_power(matrix,power) =
    (power<=1)?matrix:(power%2)?matrix * matrix_power(matrix, power-1):(let (M=matrix*matrix) matrix_power(M,power/2));

/*
  takes a transformation matrix and a list of points (a polygon) and applies the transformation matrix to the points
*/
function transform_shape(matrix, shape) =
    [ for (point=shape) [for (i=[0:2]) (matrix*concat(point,1))[i]]];


/*
 The big show!
 Takes 
 shape - a list of 3D points that define your shape. It will close your shape (connect the last point to the first point)
pitch - how high to go up each rotation
rotations - how many times to rotate    
 */
module helix_extrude(shape=[for(a=[-90:20:260]) [1.5*sin(a)+10,0,-1.5*cos(a)]], pitch = 6, rotations = 3.5, convexity = 20) {

    if(len(shape) < 3) {
        echo("<b>ERROR! Shape has fewer than three points!</b>");
    }

    /* checking that the base shape is planar
       Establish a plane from the first three points using the cross product of differences (and then making it a unit vector) and then check that points i, i+1 and i+2 are in the same plane. You might get a lot of errors, but it is better than none. 
       Non-Planar shapes don't prevent the algorithm from working, but lead to funny behavior when rendering.
       Note: opposite cross products are still on the same plane. A cross product of [0, 0, 0] means that you have three co-linear points, so we don't look for it, either.
    */
    basePlane = let (a=cross(shape[0] - shape[1], shape[0] - shape[2])) a/sqrt(a*a);
    for(i=[1:len(shape)-3]) {
        if(let(c = cross(shape[i]- shape[i+1], shape[i]-shape[i+2]), d=c/sqrt(c*c)) basePlane != d && basePlane != -d && c != [0,0,0]){
            echo("<b>WARNING! Shape not in a plane!</b>");
        } // end if
    } // end for

    /* 
    I haven't figured out a good way to use $fn, $fa, and $fs.
    Currently it uses $fn on a PER ROTATION basis if it is defined, otherwise it uses $fa.
    
    Here we figure out the number of steps - it must be an integer, and then the angle per rotation.
    */
    steps = rotations*($fn?ceil($fn):ceil(360 / $fa));
    angle = rotations * 360 / steps;

    /*
    This is the "base" transformation matrix.
    I used to use it and the matrix power function above, but that renders slower than changing the matrix each time.
    */
    MM_base = [[cos(angle), -sin(angle), 0, 0],
          [sin(angle), cos(angle), 0, 0],
          [0, 0, 1, pitch*angle/360],
          [0,0,0,1]];
    /*
    Calculates the points on our polyhedron. Takes our base shape - a set of 3D points - and runs it through a transformation matrix using the transformShape function above;
    */
    ph_points = [for(a=[0:1:steps]) for(b=[0:len(shape)-1]) let (MM = 
        [[cos(angle*a), -sin(angle*a), 0, 0],
         [sin(angle*a), cos(angle*a), 0, 0],
         [0, 0, 1, a*pitch*angle/360],
         [0, 0, 0, 1]]) transform_shape(MM,shape)[b]];
    
    /*
     List the faces.
     The ends are done first. It is assumed your shape was planar. If it wasn't, I don't know what will happen.
     Then it makes triangles:
     (a) current point, next point in this instance, current point in next instance
     (b) current point, next point in next instance, current point in next instance
    */
    ph_faces = concat([[for(a=[1:len(shape)]) len(shape)-a]], 
        [[for(a=[0:len(shape)-1]) steps*len(shape) + a]],
        [
        for(a=[0:steps-1]) 
        for(b=[0:len(shape)-1]) 
            [
                (a)*len(shape)+b, 
                (a)*len(shape)+(b+1)%len(shape),
                (a+1)*len(shape)+(b+1)%len(shape)
            ]], 
        [
        for(a=[0:steps-1]) 
        for(b=[0:len(shape)-1]) 
            [
                (a)*len(shape)+b, 
                (a+1)*len(shape)+(b+1)%len(shape), 
                (a+1)*len(shape)+b
            ]
        ]);
    /*
     And draw the polyheadron. Nothing funny here
     */
    polyhedron(points=ph_points, faces=ph_faces, convexity=convexity);
};

/* The following functions translate, scale and rotate your shape
   You must give them a vector and the shape.
   They are for your convenience.
   You can also make your own transformation matrix and use transform_shape yourself!
*/
function he_translate(tv, shape) =
    let (M = [[1, 0, 0, tv[0] ],
              [0, 1, 0, tv[1] ],
              [0, 0, 1, tv[2] ],
              [0, 0, 0, 1]]) transform_shape(M, shape);

function he_scale(sv, shape) =
    let (M = [[sv[0], 0, 0, 0 ],
              [0, sv[1], 0, 0 ],
              [0, 0, sv[2], 0 ],
              [0, 0, 0, 1]]) transform_shape(M, shape);

function he_rotate(rv, shape) =
    let (
         Mx = [[1, 0, 0, 0 ],
               [0, cos(rv[0]), -sin(rv[0]), 0 ],
               [0, sin(rv[0]), cos(rv[0]), 0 ],
               [0, 0, 0, 1]],
         My = [[cos(rv[1]), 0, -sin(rv[1]), 0 ],
               [0, 1, 0, 0 ],
               [sin(rv[1]), 0, cos(rv[1]), 0 ],
               [0, 0, 0, 1]],
         Mz = [[cos(rv[2]), -sin(rv[2]), 0, 0 ],
               [sin(rv[2]), cos(rv[2]), 0, 0 ],
               [0, 0, 1, 0 ],
               [0, 0, 0, 1]]          
         ) transform_shape(Mz*My*Mx, shape);

/*
Makes a circle of diameter d centered on the origin.
 */
function he_circle(d) =
 let(steps = $fn?$fn:max(5,min(360/$fa, PI*d/$fs)))
 [for(a=[0:steps]) [d/2 * cos(a*360/steps), d/2 * sin(a*360/steps),0]];
/*
Makes a square. Can be centered on the origin or have a vertex on the origin using center. You must give it a vector.
*/ 
function he_square(v, center=false) =
 center?[[-v[0]/2,-v[1]/2,0], [v[0]/2,-v[1]/2,0], [v[0]/2, v[1]/2, 0], [-v[0]/2, v[1]/2, 0]]:[[0,0,0], [v[0],0,0], [v[0], v[1], 0], [0, v[1], 0]];


/* START OF SECTION FOR MAKING THREADS! */

/*
 Note: for internal threads I gave another 1/8 of the pitch.
 I think this will make it work OK on 3D printers.
 The profile leaves off the inner section to avoid touching vertices and manifoldness problems. This is added back by sticking a cylinder in the middle in the fuction iso_bolt().
 Internal threads are a bit higher than the ISO standard.
 External threads are a bit deepar than the ISO standard.
 Both of these are done to prevent locking up.
 */
function iso_thread_profile(pitch = 1, diameter = 10, length=40, internal = false) =
  internal?iso_thread_profile_internal(pitch=pitch, diameter=diameter, length=length):iso_thread_profile_external(pitch=pitch, diameter=diameter, length=length);
 
 function iso_thread_profile_internal(pitch = 1, diameter = 10) =
    let(h = pitch * cos(30),
      dmin = diameter - (5/8-internal_fudge) * h,
      dmaj = diameter + (1/16+internal_fudge) * h,
      manifold_fudge_x = 3/4 * dmin/2)
     [[manifold_fudge_x,0,pitch*(1/2 + 1/16 - 1/32)],
      [dmin/2, 0, pitch*(1/2 + 1/16 - 1/32)],
      [dmaj/2, 0, pitch*15/16],
      [dmaj/2, 0, pitch],
      [dmin/2, 0, pitch*(1+1/2 - 1/16 - 1/32)],
      [manifold_fudge_x,0,pitch*(1+15/16)]];
    
 function iso_thread_profile_external(pitch = 1, diameter = 10) =
    let(h = pitch * cos(30),
      dmin = diameter - 3/4 * h,
      dmaj = diameter,
      manifold_fudge_x = 3/4 * dmin/2)
     [[manifold_fudge_x,0,0],
      [dmin/2, 0, 0],
      [dmaj/2, 0, pitch*(1/2 - 1/16 - 1/16)],
      [dmaj/2, 0, pitch*(1/2 - 1/16 + 1/16)],
      [dmin/2, 0, pitch*7/8],
      [manifold_fudge_x,0,pitch*7/8]];


/* 
THIS IS WHY YOU CAME HERE!
*/
module iso_bolt(diameter=10, pitch=1, length=40, internal=false) {
    shape = iso_thread_profile(pitch=pitch, diameter=diameter, internal=internal, length=length);
    /*
    This is to make the segments on the cylinder and the segments on the thread line up.
    It makes sure that $fa and $fs are functionally equal.
    */
    $fa = (PI*diameter/$fs < 360/$fa)?($fs*360/(PI*diameter)):$fa;
    $fs = (PI*diameter/$fs < 360/$fa)?$fs:(360/($fa*PI*diameter));
    
    cylinderD = internal?(diameter - (5/8-internal_fudge)*(pitch*cos(30))):(diameter - 3/4 * pitch*cos(30));
    difference() {
        union() {
            // we still want to shrink $fs here to the size of the bolt cylinder diameter.
            cylinder(d=cylinderD,h=length, $fs = $fs * cylinderD/diameter);
            translate([0,0,-pitch])
            helix_extrude(shape = shape, pitch=pitch, rotations=length/pitch+2);
        }
        translate([-(diameter+1)/2, -(diameter+1)/2, -4*pitch])
        cube([diameter+1, diameter+1, 4*pitch]);
        translate([-(diameter+1)/2,-(diameter+1)/2,length])
        cube([diameter+1, diameter+1, 4*pitch]);
    }
}



// a square 3x3
square = he_rotate([90,0,0], he_translate([10,0,0], he_square([3,3], center=true)));
// a circle d=3
circle = he_rotate([90,0,0], he_translate([10,0,0], he_circle(3)));

// a circle helix
echo("Circle Helix");
translate([50,0,0])
helix_extrude(circle);
// a square helix
echo("Square Helix");
translate([25,0,0])
helix_extrude(square);
// a bolt
echo("External Bolt");
iso_bolt(pitch=1,internal=false);
// a bolt with internal threading -- for subtracting!
echo("Internal Bolt");
translate([-20,0,0])
iso_bolt(pitch=1,internal=true);
