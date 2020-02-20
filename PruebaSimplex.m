%EP2 Lectura de un PPL en formato MPS
% Integrantes:
%Ocampo Giles Karen Lizeth 
%Pérez Ramos Diana Areli
%Ríos Contreras Israel
%Investigacion de Operaciones
%8B Informática

file = fopen('ejemplo1.mps','r');% Abrir el archivo y almacenarlo en la variable "file"
if file > 0% Si se abrio el archivo de manera correcta, "file" es mayor a 0
    palabra = fscanf(file, '%s', 1);% Se lee la primera palabra y se almacena en "palabra"
    if strcmp(palabra, 'NAME')%Comparamos si la palabra anterior coincide con "NAME"
        palabra = fscanf(file, '%s', 1);%Si es asi hacemos lectura de la siguiente palabra
        fprintf ('El titulo del archivo es %s \n', palabra);% Mostramos el nombre del problema
        fprintf ('NAME %s', palabra);%Mostramos lo que hemos leido hasta ahora
        palabra = fscanf(file, '%s', 1);%Hacemos otra lectura para asegurarnos que lo que sigue en el archivo es "ROWS"
        if strcmp(palabra, 'ROWS')%Si la palabra leida es "ROWS" seguimos con la lectura del archivo
            rows = [];%Creamos un vector para almacenar el nombre de los renglones (obj, r01)
            etiquetas = [];%Creamos un vector para almacenar las etiquetas (N,L)   
            while ~strcmp(palabra, 'COLUMNS')%Mientras la lectura de la palabra sea diferente de 'COLUMNS' nos movemos en el archivo y almacenamos en las matrices creadas
                palabra = fscanf(file, '%s', 1);%Volvemos a leer
                if strcmp(palabra, 'N') || strcmp(palabra, 'G') || strcmp(palabra, 'L') || strcmp(palabra, 'E') || ~strcmp(palabra, 'COLUMNS')%Si la palabra leida es igual a cualquiera de estas palabras reservadas en mps entra
                    etiquetas = [etiquetas; palabra];%Se concatena la palabra a la mattriz correspondiente
                elseif strcmp(palabra, 'COLUMNS')%De otro modo si la palabra leida es 'COLUMNS' se sale del while
                    break;
                else
                    disp('Tu archivo tiene un error en el nombre de los renglones');%De otro modo el archivo tiene un error
                    return;
                end
                palabra = fscanf(file, '%s', 1);%Volvemos a leer palabra para meter en rows
              
                if ~strcmp(palabra, 'N') || ~strcmp(palabra1, 'G') || ~strcmp(palabra1, 'L') || ~strcmp(palabra1, 'E') || ~strcmp(palabra1, 'NAME') || ~strcmp(palabra1, 'ROWS') ||  ~strcmp(palabra1, 'RHS') || ~strcmp(palabra, 'ENDATA')%Si la palabra leida es diferente de cualquier palabra reservado de la extension mps, entra   
                     rows = [rows; palabra];%Se hace lo mismo que arriba, lo que ya tenia rows se le concatena mas la palabra leida
                elseif strcmp(palabra, 'COLUMNS')%Si ya se leyo la palabra "COLUMNS" se sale del while
                    break;
                else
                    disp('Tu archivo tiene un error en las etiquetas');%Si no se tiene un error en el archivo
                end
            end
            fprintf('\nEtiquetas\n');
            disp(etiquetas);
            fprintf('\n');
            fprintf('ROWS\n');
            disp(rows);
            if strmatch('N', etiquetas)>0
                    %Se inicia el procesamiento de la primera variable
                    palabra = fscanf(file, '%s', 1); 
                    [m,n]= size(etiquetas);%Se necesita declarar una matriz de mxn(columna de la matriz A)
                    %Para poder saber esto se obtiene el tamaño de rows o de
                    %etiquetas
                    %Se declara una matriz A columna de tamaño m x 1 la llenaremos
                    %de 0
                    a = zeros(m,1);
                    %Este vector nos servira para ir creando la columna de A
                    %correspondiente a la primera columna, es decir variable x0
                    palabra = fscanf(file, '%s', 1);%Esta palabra corresponde al nombre del renglon restriccion
                    %donde se debe colocar el valor del coeficiente de la variable
                    %Por lo que se debe saber en que posicicon del vector se debera
                    %colocar para determinar el valor se vusca el nombre en el
                    %vector rows
                    A = []%Se crea la matriz A
                    i=0;%Esta i nos ayudara en los while siguientes
                    while(~strcmp(palabra,'RHS') || ~strcmp(palabra,'rhs'))%mientras el recorrido con palabra sea diferente de "RHS"
                        renglon = strmatch(palabra, rows);
                        palabra = fscanf(file, '%s', 1);
                        while(~strcmp(palabra,'obj') || ~strcmp(palabra,'rhs') )%Hasta que palabra sea igual a obj o rhs se seguira en el ciclo
                            if strcmp(palabra,'rhs')
                               break;
                            else
                            end
                            if i == 0%esta i es para la primera vez que entra pues se tienen diferentes condiciones
                                if strcmp(palabra,'RHS') || strcmp(palabra,'rhs')
                                    break;
                                elseif strcmp(palabra,'obj') || strcmp(palabra,'rhs')
                                    break;
                                else
                                    numero = str2num(palabra);
                                    a(renglon,1) = numero;
                                    i= i+1;%Aumentamos i para que entre al else desde la segunda vez en adelante
                                end
                            else
                               palabra = fscanf(file, '%s', 1);
                               if strcmp(palabra,'obj') || strcmp(palabra,'RHS') || strcmp(palabra,'rhs') 
                                   break;
                               else
                                palabra = fscanf(file, '%s', 1);
                               end
                               if strcmp(palabra,'obj')|| strcmp(palabra,'RHS') || strcmp(palabra,'rhs')
                                   break;
                               else
                                   renglon = strmatch(palabra, rows);
                                   palabra = fscanf(file, '%s', 1);
                                   if strcmp(palabra,'RHS') || strcmp(palabra,'rhs')
                                    break;
                                   elseif strcmp(palabra,'obj')|| strcmp(palabra,'rhs')
                                    break;
                                   else
                                    numero = str2num(palabra);
                                    a(renglon,1) = numero;
                                   end
                                   if strcmp(palabra,'rhs')
                                       break;
                                   else
                                       continue;
                                   end
                               end
                            end
                           
                        end
                           i=0;
                        if strcmp(palabra,'rhs') || strcmp(palabra,'RHS')
                               break;
                        else
                        end 
                        A = [A,a];
                        a = zeros(m,1);
                       
                    end
               A = [A,a];
               disp('Vector C(primera fila) y matriz A juntas'); 
               disp(A);  
               c=A(1,:);
               disp('Vector C'); 
               disp(c);
               disp('Matriz A'); 
               A = A(2:m,:)
                if (strcmp(palabra,'RHS'))
                  palabra = fscanf(file, '%s', 1);
                  columna=[]; % se crea un vector para almacenar las col
                  b=[]; % se crea un vector para almacenar los valores de b
                  while (strcmp(palabra,'rhs')) 
                  palabra = fscanf(file, '%s', 1); %itera a la sig palabra
                  columna=[columna; palabra]; % almacena la palabra en el vector columna
                  palabra = fscanf(file, '%s', 1);
                  b=[b; str2num(palabra)]; % se agrega la palabra al vector b
                  palabra = fscanf(file, '%s', 1);
                  end 
                  fprintf('Vector b\n'); 
                  disp(b); %Se muestra el vector b
                else
                    fprintf('\nTu archivo tiene un error, no contiene la palabra RHS\n');
                end 
             if (strcmp(palabra,'ENDATA'))
                disp(palabra); % se imprime la palabra "ENDATA"
                fprintf('\nSe cerrara el archivo\n');
                fclose(file); %Cerrar el archivo
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AQUI TERMINA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LA LECTURA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DEL ARCHIVO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
               %%%%%%%%%%%%%%%%%%%%%%%Empieza la parte del simplex
               %(EP2)%%%%%%%%%%%%%%%%%%%%%%%%%%
                i = 0;
                for j=2:1:size(etiquetas)%Verificamos si todas las condicinales son L, es decir menor igual que
                   if (~strcmp(etiquetas(j), 'L'))
                       i = 1;%Para esto nos sirve la i de la linea 153, la usamos como bandera, si hay algo diferente L i es igual a 1
                   end
                end
             if i == 1
                 fprintf('\nEste archivo no contiene desigualdad menor que igual que, tiene un error, vuelve a modelar');
                 return;
             end
             t = size(A);%Necesitamos saber el tamaño de la matriz A para eso es t, el tamaño
             t = t(1,1);
             %Funcion para la holgura y concatenacion%
             [A,c]=conca(A,c,t);
             XB = [];%BASICAS%
             XN = [];%NO BASICAS%
             for p=1:numel(c)%numel(x) devuelve el número de elementos, n, en array A, equivalente a prod(size(A)).(documentacion matlab)
                 if(c(p)~=0)
                     XN=[XN,p];
                 else
                     XB=[XB,p];
                 end
             end
             [zpr,bpr,XB] = metodoS(A,b,c,XB,XN);%SE manda llamar al metodo simplex
             else
                 fprintf('\nTu archivo tiene un error, no contiene la palabra ENDATA\n');
                 return;
             end
            else
                fprintf('\nTu archivo no tiene funcion objetivo\n');
                return;
            end
        else
           fprintf('\nTu archivo tiene un error, no contiene la palabra ROWS\n');
           return;
        end 
    else
        fprintf('\nTu archivo tiene un error, no contiene la palabra NAME\n');
        return;       
    end
else
    fprintf('\nNo se ha abierto el archivo\n');
    return;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FUNCIONES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%En una funcion primero se ponen las salidas de la funcion, seguido del
%nombre de la funcion y las entradas de la misma%
function [A,c] =conca(A,c,t)
%El primer paso es concatenar nuestro vector c%
ro=zeros(1,t);%Creamos un vector de ceros con tamaño de columntas t
c=[c,ro];%Las anexamos a un lado de c
h = [];
%Despues creamos la matriz identidad%
for k=1:t
    ro = zeros(1,t);
    ro(k) = 1; %conforme vamos creando la fila identidad
    h = [h;ro];%se la anexamos a h
end
   A =  [A,h];%Para que al final se la anexemos a la matriz A
end

%Esta funcion en si es el metdo simplex
function [zpr,bpr,XB] = metodoS(A,b,c,XB,XN)%funcion que se encarga de hacer el metodo simplex%
    apunt = 0;
    while 1
        apunt = apunt+1;
        fprintf('\nTabla\n');
        fprintf('Las basicas\n');
        disp(XB);%Mostramos las basicas
        fprintf('No basicas\n');
        disp(XN);%Mostramos las no basicas
        B=A(:,XB);%Como A ya tiene matriz identidad anexada, le pasamos a B todo lo que esta en las posiciones que nos indica XB
        cB=c(:,XB);%Lo mismo de arriba pero con c
        Binve=inv(B); %Sacamos la inversa de B
        cpr=c-cB*Binve*A;
        Apr=Binve*A;
        zpr=(-1)*cB*Binve*b;
        bpr=Binve*b;
        m=[cpr,zpr;Apr,bpr];
        disp(m);
        %La parte de la optimalidad%
        if(min(cpr)<0)%la funcion min nos devuelve el valo minimo de cpr
            for i=1:numel(cpr)
                if(cpr(i)==min(cpr))
                    xe=i;
                    break;
                end
            end
            acot=Apr(:,xe);%se le agregan valores a "acot"
            NA=acota(acot);%se manda a llamar la funcion de acotacion
            if(NA==numel(acot))
                disp('Problema no Acotado');
                break;
            end
            coci=bpr./Apr(:,xe);%sacamos el cociente que es coci con el operador ./ que es operacion de division de dereche de elementos
            minimo=fact(coci);%Esta es la parte de factibilidad
            for i=1:numel(coci)%Iteramos hasta el numero de elementos que tenga coci
                if(coci(i)==minimo)
                    xs=i;
                    break;
                end
            end             
            xe=XN(xe);
            xs=XB(xs);
            %Cambio de valores%
            for i=1:numel(XB)
                if(XB(i)==xs)
                    XB(i)=xe;
                    break;
                end
            end
            for i=1:numel(XN)%numel es para saber el numero de elementos en un array
                if(XN(i)==xe)
                    XN(i)=xs;
                    break;
                end
            end   
        else
            disp('Tablero Final');
            disp('XB:');
            disp(XB);
            disp('XN:');
            disp(XN);
            B=A(:,XB);
            cB=c(:,XB); 
            Binv=inv(B); 
            cpr=c-cB*Binv*A;
            Apr=Binv*A;
            zpr=-1*cB*Binv*b;
            bpr=Binv*b;
            m=[cpr,zpr;Apr,bpr];
            disp(m);
            disp('Solucion encontrada');
            break;
        end
        
    end
    fprintf('Se hicieron un total de %d iteraciones:\n',apunt);
end
 
function xs=fact(coci)
    positivos=[];
    %Crear Vector con puro positivos
    for i=1:numel(coci)
        if (coci(i)>=0)
          positivos=[positivos,coci(i)]; 
        end
    end
    xs=min(positivos);
end

function NA=acota(acot)
    for i=1:numel(acot)
        if (acot(i)>=0)
          NA=0; 
        else
          NA=i;     
        end
    end
end
    