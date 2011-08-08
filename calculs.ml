module Propu = struct
      let table = 
	[(0.,0.); (0.05, 893.); (0.5, 798.); (1.,739.); 
	 (1.5,658.); (2.,586.); (2.5, 513.); (2.97, 417.); 
	 (3.2, 225.); (3.47, 67.); (3.59, 0.)] 
      let propu t = 	
	let rec find (t0, p0) = function
	  | [] -> 0.
	  | (t1, p1)::r when t > t1 -> find (t1, p1) r
	  | (t1, p1)::r  ->
	    let a = (p1-.p0)/.(t1-.t0) in
	    let b = p0 -.a*.t0 in
	    a*.t+.b
	in
	find (List.hd table) (List.tl table)
end
module Masse = struct
     let tf = 3.59
     let m0, mf = 1.685, 0.652
     let masse_vide = 6.55
     let m_propu t =   
       if t > tf then 0. 
       else
	 let a = (mf-.m0)/.tf in
	 a*.t+.m0

     let centre_masse = 1.36
     let masse t = masse_vide +. m_propu t
     let masse' t = if t > tf then 0. else (mf-.m0)/.tf
end
module Const = struct
  let cx = 0.7
  let rho_air = 1.2
  let pi = 3.141
  let d = 0.08 
  let s_aileron = 0.000001 *. 150.*.12.
  let s = (pi*.d*.d)/.4. +. s_aileron*.4.
  let k = 0.5*.cx*.rho_air*.s
    
  let g = 9.81
  let pas = 0.00001
end;;


let rad a = Const.pi *. a /. 180.;;
let longueur_rampe = 4.;;

let l0 = 0.2;;
let masse_ressort = 0.5;;

module Data = struct
  let theta = ref (rad 80.)
  let theta' = ref 0.
  let theta'' = ref 0.

  let r = ref Masse.centre_masse
  let r' = ref 0.
  let r'' = ref 0.

  let xe = ref (l0/.2.)
  let xe' = ref 0. 
  let xe'' = ref 0.

  let t = ref 0.

  let init () =
    theta := rad 80.; theta' := 0.; theta'' := 0.;
    r :=  Masse.centre_masse; r' := 0.; r'' := 0.;
    xe := l0/.2.; xe' := 0.; xe'' := 0.;
    t:= 0.

  let update ntheta ntheta' ntheta'' nr nr' nr'' nt = 
    theta := ntheta; theta' := ntheta'; theta'' := ntheta'';
    r := nr; r' := nr'; r'' := nr'';
    t := nt;
  ;;
  let update_ressort nx nx' nx'' = 
    xe := nx; xe' := nx'; xe'' := nx''


  let h () = !r*.sin !theta
  let v () = sqrt (!r'**2. +. (!r*.(!theta'))**2.)
  let acc () = 
    sqrt ((!r''-.(!r')*.(!theta')**2.)**2. +. (2.*.(!r')*.(!theta')+.(!r')*.(!theta''))**2.)
  let acc_r () = 
    !r'' -. !r' *. (!theta')**2.

  let print () = 
    Printf.printf "t = %f \t h = %f v = %f \t a = %f \t a_r = %f \nr = %f \t r' = %f\tr'' = %f\ntheta = %d° \t theta' = %f theta'' = %f \n\n" 
      !t (h ()) (v ()) (acc ()) (acc_r ()) !r !r' !r'' (int_of_float ((!theta*.180.)/.Const.pi)) !theta' !theta'';;

  let print_ressort() = 
    Printf.printf "t = %f; x = %f; x' = %f; x'' = %f\n" !t !xe !xe' !xe''
end
;;

let k = 1.22501;;
let rampe () = 
  while !Data.r -.Masse.centre_masse < longueur_rampe do
    if !Data.xe < l0/.2. || !Data.xe > l0 then failwith "mauvais k"
    else 
      begin
	let nt = !Data.t +. Const.pas in
	let r'' = (Propu.propu nt -. Const.k*.(!Data.r')**2. -. !Data.r'*.Masse.masse' nt)/.(Masse.masse nt) 
	  -. Const.g*.sin (!Data.theta) in
	let r' = !Data.r' +. r'' *. Const.pas in
	let r = !Data.r +. r' *. Const.pas +. 0.5*.r''*.(Const.pas**2.)  in

	Data.update !Data.theta 0. 0. r r' r'' nt;

	let x'' = (Propu.propu nt -. k*.(2.*.(!Data.xe)-.l0))/.masse_ressort -. Data.acc_r () -. Const.g in
	let x' = !Data.xe' +. x'' *. Const.pas in
	let x = !Data.xe +. x' *. Const.pas +. 0.5*.x''*.(Const.pas**2.)  in
    
	Data.update_ressort x'' x' x;
	Data.print_ressort();
      end
  done;
  print_endline "Sortie de la rampe";
;;

let poussee () = 
  let old_h = ref 0. in
  while !old_h < Data.h () do
 
	old_h := Data.h ();
	let nt = !Data.t +. Const.pas in
	let theta'' = (-.Const.g *. cos (!Data.theta) -. 2.*.(!Data.r')*.(!Data.theta'))/.(!Data.r) in
	let r_air = Const.k *. (!Data.r'**2.+. (!Data.r*.(!Data.theta'))**2.) in
	let r'' = !Data.r*.(!Data.theta'**2.) 
	  +. (Propu.propu nt -. r_air)/.(Masse.masse nt)
	  -. Const.g*.sin(!Data.theta)
	in
    
	let theta' = !Data.theta' +. theta'' *. Const.pas in
	let theta = !Data.theta +. theta' *. Const.pas +. 0.5*.theta''*.(Const.pas**2.)  in
	
	let r' = !Data.r' +. r'' *. Const.pas in
	let r = !Data.r +. r' *. Const.pas +. 0.5*.r''*.(Const.pas**2.)  in
 
	Data.update theta theta' theta'' r r' r'' nt;
 
  done;
  print_endline "apogée"
;;


  Data.init();
  rampe ();;
  poussee ();
Data.print();;
