type propu = {nom:string; poussee:(float*float) list; m:float; m_vide:float};; 
module Propu = struct
      let tables = 
	[{nom="cariacou"; poussee= [(0., 0.); (0.02, 320.); (0.04, 170.); (0.06, 217.); (0.94, 85.); (0.95, 0.)];
	  m=0.220; m_vide=0.150};
	 {nom="barasinga"; poussee= [(0.,0.); (0.05, 893.); (0.5, 798.); (1.,739.); (1.5,658.); (2.,586.); (2.5, 513.); (2.97, 417.); (3.2, 225.); (3.47, 67.); (3.59, 0.)]; m=1.685; m_vide=0.652}
	 ]
      let p = ref "cariacou"
      let car () = List.find (fun e -> e.nom = !p) tables
      let propu t = 
	let table = (car ()).poussee in
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
     let masse_vide = ref 6.55
     let m_propu t =   
       let m0, mf = (Propu.car()).m,  (Propu.car()).m_vide in 
       if t > tf then 0. 
       else
	 let a = (mf-.m0)/.tf in
	 a*.t+.m0

     let centre_masse = 1.36
     let masse t = !masse_vide +. m_propu t
end


module Const = struct
  let cx = ref 0.45
  let rho_air h  = 1.225 *. ((20000.-.h)/.(20000.+.h))
  let pi = 3.141
  let d = 0.08 
  let s_aileron = 0.000001 *. 150.*.12.
  let s = ref ((pi*.d*.d)/.4. +. s_aileron*.4.)
  let k h () = 0.5*.(!cx)*.(rho_air h)*.(!s)
    
  let g = 9.81
  let pas = 0.001

  let longueur_rampe = ref 4.;;
  let l_potar = ref 0.08
  let l0 = ref 0.2
  let masse_ressort = ref 0.6
  let k_ressort = ref 14000.;;

end


let rad a = Const.pi *. a /. 180.;;
let debug = ref false;;
module Data = struct
  let theta = ref (rad 80.)
  let x = ref 0.
  let x' = ref 0.
  let x'' = ref 0.

  let z = ref 0.
  let z' = ref 0. 
  let z'' = ref 0.

  let t = ref 0.


  let xr = ref (!Const.l0/.2.)
  let xr' = ref 0.
  let xr'' = ref 0.

  let init () =
    theta := rad 80.;
    x :=  0.; x' := 0.; x'' := 0.;
    z := 0.; z' := 0.; z'' := 0.;
    xr := !Const.l0/.2.; xr' := 0.; xr'' := 0.;
    t:= 0.

  let update ntheta nx nx' nx'' nz nz' nz'' nt = 
    theta := ntheta; 
    x := nx; x' := nx'; x'' := nx'';
    z := nz; z' := nz'; z'' := nz'';
    t := nt;
  ;;

  let update_ressort nxr nxr' nxr'' =
    xr := nxr; xr' := nxr'; xr'' := nxr''
  ;;

  let v () = sqrt (!z'**2. +. !x' **2.)
  let a () = !x'' *. cos !theta +. !z'' *. sin !theta
  let l () = 2.*.(!Const.l0-.(!xr));;
  let print () = 
    Printf.printf "t : %f \t h = %f \t v = %f \t theta = %d° \t a = %f \n" !t (!z) (v())  (int_of_float ((!theta*.180.)/.Const.pi)) (a () )

  let print_ressort () = 
    Printf.printf "x = %f, x'' = %f, l = %f, F = %f\n" !xr !xr'' (l()) (!Const.k_ressort*.(2.*.(!xr)-.(!Const.l0)))

 
end


exception Potar;;
let rampe () = 
  while !Data.z < !Const.longueur_rampe *. sin (!Data.theta) do
    if abs_float (Data.l() -. !Const.l0) > !Const.l_potar then raise Potar;
    let nt = !Data.t +. Const.pas in
    let a = !Data.theta in
    let m = Masse.masse nt in
    let h = !Data.z in
    let x'' = (cos a /. m) *. (Propu.propu nt -. (Const.k h ())*.Data.v()**2. -. m*.Const.g*.sin a) in

    let z'' = (sin a /. m) *. (Propu.propu nt -. (Const.k h ())*.Data.v()**2. -. m*.Const.g) *.sin a in

    let x' = !Data.x' +. x'' *. Const.pas in
    let x = !Data.x +. x' *. Const.pas +. 0.5*.x''*.(Const.pas**2.)  in

    let z' = !Data.z' +. z'' *. Const.pas in
    let z = !Data.z +. z' *. Const.pas +. 0.5*.z''*.(Const.pas**2.)  in   
    
    let xr'' = ((Propu.propu nt +. !Const.k_ressort *. (!Const.l0 -. 2.*.(!Data.xr)))/.(!Const.masse_ressort)) -. Const.g *. sin a -. Data.a () in
    let xr' = !Data.xr' +. xr'' *. Const.pas in
    let xr = !Data.xr +. xr' *. Const.pas +. 0.5*.xr''*.(Const.pas**2.)  in

    Data.update a x x' x'' z z' z'' nt;
    Data.update_ressort xr xr' xr'';
    if !debug then
      begin 
	Data.print();
	Data.print_ressort();
	print_newline();
      end
  done;
  if !debug then print_endline "Sortie de la rampe";
;;

let poussee () = 
  while !Data.z' > 0. do
    if abs_float (Data.l() -. !Const.l0) > !Const.l_potar then raise Potar;

    let nt = !Data.t +. Const.pas in
    
    let a = !Data.theta in
    let m = Masse.masse nt in
    let h = !Data.z in

    let x'' = (cos a /. m) *. (Propu.propu nt -. (Const.k h ())*.Data.v()**2.) in
    let z'' = (sin a /. m) *. (Propu.propu nt -. (Const.k h ())*.Data.v()**2.) -. Const.g in

    let x' = !Data.x' +. x'' *. Const.pas in
    let x = !Data.x +. x' *. Const.pas +. 0.5*.x''*.(Const.pas**2.)  in

    let z' = !Data.z' +. z'' *. Const.pas in
    let z = !Data.z +. z' *. Const.pas +. 0.5*.z''*.(Const.pas**2.)  in   
    

    let xr'' = ((Propu.propu nt +. !Const.k_ressort *. (!Const.l0 -. 2.*.(!Data.xr)))/.(!Const.masse_ressort)) -. Const.g *. sin a -. Data.a () in
    let xr' = !Data.xr' +. xr'' *. Const.pas in
    let xr = !Data.xr +. xr' *. Const.pas +. 0.5*.xr''*.(Const.pas**2.)  in


    let theta = atan (z'/.x') in

    
    Data.update theta x x' x'' z z' z'' nt;
    Data.update_ressort xr xr' xr'';
    if !debug then 
      begin 
	Data.print();
	Data.print_ressort(); 
	print_newline(); 
      end;
  done;
  if !debug then print_endline "apogée"; 
;;

let descente () = 
  while !Data.z > 0. do
    if abs_float (Data.l() -. !Const.l0) > !Const.l_potar then raise Potar;

    let nt = !Data.t +. Const.pas in
    
    let a = !Data.theta in
    let m = Masse.masse nt in
    let h = !Data.z in

    let x'' = (cos a /. m) *. (Propu.propu nt -. (Const.k h ())*.Data.v()**2.) in
    let z'' = (sin a /. m) *. (Propu.propu nt +. (Const.k h ())*.Data.v()**2.) -. Const.g in

    let x' = !Data.x' +. x'' *. Const.pas in
    let x = !Data.x +. x' *. Const.pas +. 0.5*.x''*.(Const.pas**2.)  in

    let z' = !Data.z' +. z'' *. Const.pas in
    let z = !Data.z +. z' *. Const.pas +. 0.5*.z''*.(Const.pas**2.)  in   
    

    let xr'' = ((Propu.propu nt +. !Const.k_ressort *. (!Const.l0 -. 2.*.(!Data.xr)))/.(!Const.masse_ressort)) -. Const.g *. sin a -. Data.a () in
    let xr' = !Data.xr' +. xr'' *. Const.pas in
    let xr = !Data.xr +. xr' *. Const.pas +. 0.5*.xr''*.(Const.pas**2.)  in


    let theta = atan (z'/.x') in

    
    Data.update theta x x' x'' z z' z'' nt;
    Data.update_ressort xr xr' xr'';
    if !debug then 
      begin 
	Data.print();
	Data.print_ressort(); 
	print_newline(); 
      end;
  done;
  if !debug then print_endline "impact"; 
;;


let rec search (a, b) = 
  let m = (a+.b)/.2. in
  if b-.a < 1. then b
  else
    try
      Const.k_ressort := m;
      Data.init();
      rampe();
      poussee();
      descente();
      search (a, m)
    with Potar -> 
      search(m, b)
;;

let () = 
  Arg.parse 
    [("-mvide", Arg.Float(fun m -> Masse.masse_vide := m), "masse à vide de la fusée");
     ("-mressort", Arg.Float(fun m -> Const.masse_ressort := m), "masse du ressort");
     ("-l0", Arg.Float(fun l -> Const.l0 := l), "masse du ressort");
     ("-cx", Arg.Float(fun cx -> Const.cx := cx), "Cx de la fusee");
     ("-l_potar", Arg.Float(fun l -> Const.l_potar := l), "longueur du potentiometre");
     ("-propu", Arg.String(fun n -> Propu.p := n), "nom du propu : [barasinga, cariacou]");
     ("-s", Arg.Float(fun s -> Const.s := s), "surface au vent");
    ] 
    (fun s -> ()) "Programme de calcul de Mercury";  
  let k = search (1000., 400000.) in
  Const.k_ressort := k;
  Printf.printf "k = %f N.m^-1\n" !Const.k_ressort;
  (* debug := true; *)
  Data.init();
  rampe();
  Data.print();
  poussee();
  Data.print();
  descente();
  Data.print();;


