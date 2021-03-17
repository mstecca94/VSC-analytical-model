function [ Semicond ] = semiconductor_data ( Semicond  )
%% here I can add different switches and switch to compare them
topology = Semicond.topology ;

Type = Semicond.Type ;

switch topology
    case '2LC'
        switch Type
            case 'SKiiP_38GB12E4V1'
            %% IGBT parameters

            Semicond.Param.U__Bt = 600;

            Semicond.Param.rt = 0.0046;
            Semicond.Param.Vt = 0.5289575289575286;

            %% Diode parameters

            Semicond.Param.U__Bd = 600;

            Semicond.Param.rd = 0.0032;
            Semicond.Param.Vd = 0.5435767716801709;

            %% Switching functions 

            % E reverse recovery
            err=[590.8997622579401, 0.02731178753581371
            519.3820738422775, 0.025940198727978384
            456.1152540995266, 0.02476520431593278
            375.43535244752405, 0.02299798833641517
            308.50181862516, 0.021724544327718275
            239.80960315363836, 0.017959990246479582
            185.7598601995408, 0.015386584845467656
            149.09881535366674, 0.014202751305549335
            0, 0
            ];

            % E turn off
            eoff=[582.6367016845143, 0.05751386828683479
            525.8986446670596, 0.0528483327576047
            448.1112714120253, 0.04649543819722431
            352.92300814825353, 0.03915164692256111
            266.01304533355005, 0.031107430963363332
            200.1798305325829, 0.023853554954991564
            159.92827098530879, 0.02007955215085444
            149.86690508605452, 0.01908621705646881
            0, 0
            ];


            % E turn on
            eon=[155.6031942779347, 0.01150956048198646
            238.96836201816598, 0.015468575376424917
            354.34336455814525, 0.022705992319102676
            421.1671712758825, 0.027567512649096786
            489.8258590209904, 0.03242842338406518
            554.8178326864853, 0.0371908845223823
            595.1059679353017, 0.03976886188608701
            0, 0 
            ];

            RjcT = 0.17 ;
            RjcD = 0.26 ;
            Rcs = 0.05 ;
            weight = 0.076 ;
            Semicond.Param.Rcs = Rcs;        

            case 'FF300R12KE4'
            %% IGBT parameters
            Semicond.Param.U__Bt = 600;
            Semicond.Param.rt = 0.0043;
            Semicond.Param.Vt = 0.8391608391608392;
            %% Diode parameters
            Semicond.Param.U__Bd = 600;
            Semicond.Param.rd = 0.0030;
            Semicond.Param.Vd = 0.7736254773291813;

            %% Switching functions 

            % E reverse recovery
            err=[581.3343108504398, 0.03340175953079179
            496.3812316715543, 0.03246334310850439
            395.47800586510266, 0.030527859237536654
            338.54545454545456, 0.02917888563049853
            256.7008797653959, 0.026070381231671554
            189.88269794721407, 0.022434017595307918
            130.07624633431087, 0.018563049853372433
            87.24926686217009, 0.015102639296187687
            54.419354838709694, 0.011700879765395897
            0, 0
            ];

            % E turn off
            eoff=[591.8498715186751, 0.06758052794144366
            560.1422379110752, 0.06478041373582163
            511.6152308770471, 0.0601118695979443
            439.33085887818925, 0.052991668180756346
            369.04768084719797, 0.045753497547174705
            279.98286915669524, 0.0361789134892413
            203.78955018558412, 0.0278883640044644
            138.48210345991123, 0.02076543722584162
            86.07729644145665, 0.014457782853583206
            39.65271108573215, 0.00814779245723779
            0, 0
            ];


            % E turn on
            eon=[591.6240558569314, 0.07097905365068652
            557.1754873205803, 0.06443001012277104
            497.1630285254496, 0.052617514989487885
            442.81153476782515, 0.04560749604173696
            379.5208555039322, 0.03813221896332442
            276.5333402548862, 0.028094323461468576
            173.4835310301866, 0.018993952293196975
            113.05837464635191, 0.013392555869909409
            37.76053157525892, 0.006625094089859063
            0, 0 
            ];

            RjcT = 0.093 ;
            RjcD = 0.15 ;
            
            weight = 0.340 ;
            
            RcsT = 0.032 ;
            RcsD = 0.052 ;
            
            Semicond.Param.RcsT = RcsT;        
            Semicond.Param.RcsD = RcsD;        
         case 'SK250MB120SCTE2'
            %% IGBT parameters

            Semicond.Param.U__Bt = 600;

            Semicond.Param.rt = 0.0127;
            Semicond.Param.Vt = 0;

            %% Diode parameters

            Semicond.Param.U__Bd = 600;

            Semicond.Param.rd = 0.0107;
            Semicond.Param.Vd = 0;

            %% Switching functions 

            % E reverse recovery
            err=[23.864941372676682, 0.0005513455391299906
            80.1032835689886, 0.0010539386004052367
            170.65781100277223, 0.0018764766551338258
            233.55838228086986, 0.0023324750751585466
            307.8918302111754, 0.002833678343366664
            390.7966440154779, 0.0033342232885429854
            452.72509161662185, 0.0035829231005551857
            499.4133610316652, 0.003809751958511017
            0, 0
            ];

            % E turn off
            eoff=[32.423140786037706, 0.0004124393794208269
            60.055152840663894, 0.000548565952995737
            107.74846208425073, 0.0013283130106575211
            152.58464936983853, 0.0021082795093299005
            199.33655667795568, 0.0030033062445596944
            240.3566647892269, 0.0037144414129075213
            282.3488965774517, 0.004632875189267872
            323.39972643020667, 0.005666588643196231
            372.09588109223097, 0.006976212594451069
            415.04926450687947, 0.007986738448259468
            453.2385835814234, 0.008974588730972635
            496.191966996072, 0.009985114584781032
            0, 0
            ];


            % E turn on
            eon=[22.010664833115143, 0.0010814418737336452
            123.21247010116228, 0.003700397188229187
            208.16686294445947, 0.005721522042849517
            270.2138086913269, 0.007214452384958051
            330.34942323587705, 0.008638405102734966
            388.57151216800395, 0.00997033889006737
            425.81504048686634, 0.011027386238122754
            470.69511597457404, 0.012268178859053038
            500.29770830437917, 0.013095398322007738
            0, 0 
            ];
            
            RjcT = 0 ;
            RjcD = 0 ;
            Rcs = 0.15;
            Semicond.Param.Rcs = Rcs;        

            weight = 0.035 ;        
        
            case 'SKM400GB125D'
            %% IGBT parameters

            Semicond.Param.U__Bt = 600;

            Semicond.Param.rt = 0.0072;
            Semicond.Param.Vt = 0.35;

            %% Diode parameters

            Semicond.Param.U__Bd = 600;

            Semicond.Param.rd = 0.0035;
            Semicond.Param.Vd = 0.37;

            %% Switching functions 

            % E reverse recovery
            err=[77.38957534441127, 0.007595512001136215
            145.94991241774366, 0.011210528807461075
            224.49936088623778, 0.014228092600482897
            298.5229370828006, 0.016739099559721637
            390.40382521422134, 0.019073048335937134
            475.4485631775789, 0.02022345310798656
            549.2448989253419, 0.020700658050466322
            0, 0
            ];

            % E turn off
            eoff=[591.5530071463203, 0.043218690052069075
            523.7325294601261, 0.03701309104707297
            418.03655374496714, 0.028339165630204277
            351.1647026582626, 0.024421970192044194
            268.69743027594495, 0.019483420115100035
            227.50806332609574, 0.016420937243080302
            161.80618504542866, 0.011826106204017965
            115.01359698125935, 0.008846786263886838
            76.01871956489663, 0.00637814364314776
            0, 0
            ];


            % E turn on
            eon=[78.79089144534399, 0.010137291104483276
            125.01538607205416, 0.013846518013539748
            167.87861572693274, 0.01747242342470294
            222.03758935757233, 0.02219523741892724
            283.0137764522085, 0.027932111915921043
            350.86398712304117, 0.035191497419874075
            418.8372863703072, 0.043552525682904894
            466.7376793069167, 0.05226104246555887
            0, 0 
            ];

            RjcT = 0.05 ;
            RjcD = 0.125 ;
            Rcs = 0.038;
            Semicond.Param.Rcs = Rcs;        

            weight = 0.325 ;

            case 'SKM300GB12E4'
            %% IGBT parameters

            Semicond.Param.U__Bt = 600;

            Semicond.Param.rt = 0.0051;
            Semicond.Param.Vt = 0.58;

            %% Diode parameters

            Semicond.Param.U__Bd = 600;

            Semicond.Param.rd = 0.0037;
            Semicond.Param.Vd = 0.41;

            %% Switching functions 

            % E reverse recovery
            err=[102.98694712115955, 0.012561727571203593
            148.09545418145905, 0.015379441048501183
            223.87021861570088, 0.01977898192883644
            282.47604302912214, 0.021877585884739592
            345.61033553105733, 0.02504017297066438
            410.49035047913947, 0.02571283666550994
            473.5805995248646, 0.02691989429570507
            526.7490590716174, 0.02759949817152922
            575.4130742332434, 0.028281771347729745
            0, 0
            ];

            % E turn off
            eoff=[102.13410565090891, 0.014695566292074844
            218.67309078290583, 0.029026506152737358
            260.23009369244323, 0.03415743533619836
            347.8205162426928, 0.043172196567279705
            392.07217788217713, 0.04794597336038221
            461.56474388062895, 0.05341590369164241
            512.118623708726, 0.05800816805915168
            560.8426981288204, 0.06135707231134717
            588.8303125750742, 0.06400715372500865
            0, 0
            ];


            % E turn on
            eon=[101.18116541654433, 0.012385019886287804
            142.67810906761338, 0.01484931799375383
            210.3448736086272, 0.01925366361476656
            274.396070789846, 0.023126818460881396
            324.95395456850764, 0.027896858233457
            369.1976083068628, 0.032315084216426856
            423.4471345060459, 0.040994047460160665
            456.01927234871744, 0.047196967674772425
            493.09585457651553, 0.05339721858900781
            523.886234417959, 0.060490083549101764
            554.6325708031924, 0.06562741905346606
            579.9695699757093, 0.07059018231321568
            0, 0 
            ]; 
            
            RjcT = 0.11 ;
            RjcD = 0.17 ;
            Rcs = 0.038;
            Semicond.Param.Rcs = Rcs;        

            weight = 0.325 ;
            
            case 'SEMiX303GB12Vs'
            %% IGBT parameters

            Semicond.Param.U__Bt = 600;

            Semicond.Param.rt = 0.0050;
            Semicond.Param.Vt = 0.415;

            %% Diode parameters

            Semicond.Param.U__Bd = 600;

            Semicond.Param.rd = 0.0048;
            Semicond.Param.Vd = 0.42;

            %% Switching functions 

            % E reverse recovery
            err=[99.49329079478278, 0.01239889274655158
            201.5783053392137, 0.018173031810077914
            287.6794595101811, 0.021378905883456895
            360.6624753683026, 0.023674580088204958
            457.91123205404904, 0.025078821431922693
            525.2134747114573, 0.025413812517594087
            579.4426198742611, 0.0260434456225955
            0, 0
            ];

            % E turn off
            eoff=[102.49413531012482, 0.017671483531950846
            160.66435206906257, 0.023724781833536657
            244.1475086797411, 0.03189922116918459
            333.2663976728911, 0.04082950173594822
            385.8459228675988, 0.04733180069437929
            431.8888993150042, 0.053981420662475375
            464.7743267336024, 0.05866660410997467
            508.9312189171437, 0.06486347001970537
            556.8377592192924, 0.07136342310218635
            582.2351506052361, 0.07574364267617528
            0, 0
            ];


            % E turn on
            eon=[99.44261987426106, 0.011043445622595488
            157.54527540583658, 0.015289481092239848
            230.65778361640238, 0.021049075724875688
            309.4341747208408, 0.02831753776860281
            392.0784460917706, 0.03905179694097778
            441.9780426011072, 0.04886600356573146
            492.8685371117575, 0.06018673172562636
            542.8244346438961, 0.07150699071033124
            574.8709768227457, 0.07875199399455757
            0, 0 
            ];
                    
            RjcT = 0.1 ;
            RjcD = 0.19 ;
            Rcs = 0.04 ;
            Semicond.Param.Rcs = Rcs;        

            weight = 0.3 ;
            
            case 'SKM300GB12T4'
            %% IGBT parameters

            Semicond.Param.U__Bt = 600;
            Semicond.Param.b__0t = 0; 
            Semicond.Param.b__1t = 1.7328e-04;
            Semicond.Param.b__2t = 4.8252e-08; 

            Semicond.Param.rt = 0.0050;
            Semicond.Param.Vt = 0.56;

            %% Diode parameters

            Semicond.Param.U__Bd = 600;
            Semicond.Param.b__0d = 0;
            Semicond.Param.b__1d = 0.00011336152231240664;
            Semicond.Param.b__2d = -1.1406623528819162e-7;

            Semicond.Param.rd = 0.0038;
            Semicond.Param.Vd = 0.50;

            %% Switching functions 

            % E reverse recovery
            err=[575	0.02824644549763032
            443.2692307692308	0.026729857819905195
            315.3846153846154	0.023696682464454943
            200.9615384615385	0.018388625592417038
            100.96153846153851	0.012511848341232201
            0, 0
            ];

            % E turn off
            eoff=[542.3076923076922	0.04568720379146918
            448.07692307692304	0.0396208530805687
            338.4615384615384	0.03241706161137439
            219.23076923076928	0.021990521327014192
            152.88461538461542	0.01611374407582937
            101.92307692307696	0.010995260663507075
            0, 0
            ];


            % E turn on
            eon=[562.5	0.06767772511848341
            483.65384615384613	0.05251184834123221
            372.1153846153846	0.03336492890995259
            286.53846153846155	0.02483412322274879
            153.84615384615392	0.014028436018957327
            100.96153846153851	0.010426540284360158
            0, 0 
            ];
        
            RjcT = 0.11 ;
            RjcD = 0.17 ;
            Rcs = 0.038;
            Semicond.Param.Rcs = Rcs;        

            weight = 0.325 ;

            case 'SKM200GB12T4'
            %% IGBT parameters

            Semicond.Param.U__Bt = 600;

            Semicond.Param.rt = 0.0073;
            Semicond.Param.Vt = 0.53;

            %% Diode parameters

            Semicond.Param.U__Bd = 600;

            Semicond.Param.rd = 0.0057;
            Semicond.Param.Vd = 0.43;

            %% Switching functions 

            % E reverse recovery
            err=[53.13481920453248, 0.006569037656903781
            113.75054946844008, 0.009832635983263616
            216.37985770802467, 0.013472803347280349
            280.9550168503655, 0.014100418410041851
            353.3122771600215, 0.014728033472803354
            397.67025381371786, 0.014728033472803354

            0, 0
            ];

            % E turn off
            eoff=[394.75277990329346, 0.03719665271966529
            351.26744053530444, 0.03355648535564854
            284.5009198508702, 0.027405857740585786
            220.82770297771188, 0.022008368200836834
            140.84789085521706, 0.015230125523012561
            100.47539195415393, 0.011589958158995828
            50.80344496361292, 0.006443514644351472

            0, 0
            ];


            % E turn on
            eon=[57.01610146036495, 0.00694560669456068
            106.70107288800605, 0.011589958158995828
            153.27971606727115, 0.015983263598326372
            216.10308841964732, 0.02414225941422595
            271.10854240268304, 0.03368200836820085
            302.8328150692737, 0.040711297071129714
            350.75948748840017, 0.05313807531380754
            393.23868909040596, 0.06556485355648536
            0, 0 
            ];

            RjcT = 0.14 ;
            RjcD = 0.26 ;
            Rcs = 0.038;
            Semicond.Param.Rcs = Rcs;        
            
            weight = 0.325 ;
        case 'cas300m12bm2'
            %% IGBT parameters

            Semicond.Param.U__Bt = 600;

            Semicond.Param.rt = 0.0073;
            Semicond.Param.Vt = 0.53;

            %% Diode parameters

            Semicond.Param.U__Bd = 600;

            Semicond.Param.rd = 0.0057;
            Semicond.Param.Vd = 0.43;

            %% Switching functions 

            % E reverse recovery
            err=[53.13481920453248, 0.006569037656903781
            113.75054946844008, 0.009832635983263616
            216.37985770802467, 0.013472803347280349
            280.9550168503655, 0.014100418410041851
            353.3122771600215, 0.014728033472803354
            397.67025381371786, 0.014728033472803354

            0, 0
            ];

            % E turn off
            eoff=[394.75277990329346, 0.03719665271966529
            351.26744053530444, 0.03355648535564854
            284.5009198508702, 0.027405857740585786
            220.82770297771188, 0.022008368200836834
            140.84789085521706, 0.015230125523012561
            100.47539195415393, 0.011589958158995828
            50.80344496361292, 0.006443514644351472

            0, 0
            ];


            % E turn on
            eon=[57.01610146036495, 0.00694560669456068
            106.70107288800605, 0.011589958158995828
            153.27971606727115, 0.015983263598326372
            216.10308841964732, 0.02414225941422595
            271.10854240268304, 0.03368200836820085
            302.8328150692737, 0.040711297071129714
            350.75948748840017, 0.05313807531380754
            393.23868909040596, 0.06556485355648536
            0, 0 
            ];

            RjcT = 0.07 ;
            RjcD = 0.073 ;
            Rcs = 0.038;
            
            weight = 0.3 ;
            Semicond.Param.Rcs = Rcs;        

        end
        
        Semicond.Param.eon=sort(eon);
        Semicond.Param.eoff=sort(eoff);
        Semicond.Param.err=sort(err);

        eonx= pchip(eon(:,1),eon(:,2),0:600);
        eoffx= pchip(eoff(:,1),eoff(:,2),0:600);
        errx= pchip(err(:,1),err(:,2),0:600);

        Semicond.Param.eonx=eonx;
        Semicond.Param.eoffx=eoffx;
        Semicond.Param.errx=errx;
        Semicond.Param.RjcT = RjcT;
        Semicond.Param.RjcD = RjcD;
        Semicond.Param.weight = weight ;

    case '3L_NPC'
      switch Type
          
            case 'SEMiX305MLI07E4_NPCC_650V'
                %% Diode 1  and Switch 2
                % IGBT parameters

                Semicond.Param.U__Bt2 = 300;

                Semicond.Param.rt2 = 0.0038;
                Semicond.Param.Vt2 = 0.735;

                % Diode parameters

                Semicond.Param.U__Bd1 = 300;

                Semicond.Param.rd1 = 0.0024;
                Semicond.Param.Vd1 = 0.97;

                % E reverse recovery
                err_d1=[54.463641826923094, 0.004230769230769229
                141.80438701923075, 0.005865384615384616
                256.32361778846143, 0.007403846153846153
                358.63731971153845, 0.00875
                459.0715144230769, 0.010048076923076923
                538.8762019230769, 0.011298076923076921
                0, 0
                ];

                % E turn off
                eoff_t2=[55.401141826923094, 0.004230769230769229
                128.82812500000003, 0.007451923076923075
                203.12499999999997, 0.009951923076923077
                279.29236778846155, 0.012403846153846154
                346.0712139423076, 0.014711538461538463
                413.8100961538462, 0.017259615384615387
                475.01352163461524, 0.020096153846153847
                523.0784254807692, 0.022788461538461542
                552.2896634615383, 0.024375000000000004
                0, 0
                ];


                % E turn on
                eon_t2=[51.371694711538424, 0.0012499999999999976
                136.74729567307693, 0.0019230769230769197
                241.8194110576923, 0.0026923076923076883
                339.4050480769231, 0.0036057692307692284
                453.0048076923075, 0.005336538461538462
                535.5679086538461, 0.006009615384615384
                0, 0 
                ];
            
                %% Diode 5  and Switch 1
                % IGBT parameters

                Semicond.Param.U__Bt1 = 300;

                Semicond.Param.rt1 = 0.0038;
                Semicond.Param.Vt1 = 0.735;

                % Diode parameters

                Semicond.Param.U__Bd5 = 300;

                Semicond.Param.rd5 = 0.0025;
                Semicond.Param.Vd5 = 1;

                % E reverse recovery
                err_d5=[53.58009447371825, 0.002167152435785153
                128.04060359804163, 0.0029769271633476417
                240.15994486640153, 0.004027337075909198
                348.5577682378786, 0.004795474450458723
                452.31661617539385, 0.005281195709916864
                546.8492009935534, 0.005483782968886838
                0, 0
                ];

                % E turn off
                eoff_t1=[53.17951442231765, 0.004983991155651912
                127.2652873695244, 0.00875145371792847
                191.9783485764333, 0.011437637295582132
                261.2614682192135, 0.014030639348734362
                314.87889273356404, 0.01624560295194475
                381.4699421384371, 0.019166822208502632
                433.3213685767205, 0.021757096297147116
                478.7161337564072, 0.024158566526439718
                524.1970451837069, 0.02702953380522333
                552.9828138236013, 0.028911972892647412
                0, 0
                ];

                % E turn on
                eon_t1=[55.22548780312712, 0.0011345460810636178
                149.8097603698546, 0.0016188315697282032
                249.00716449625986, 0.0022446840586368734
                372.16184008384903, 0.003437666010998
                489.82899969848813, 0.004723685910781199
                548.6496575686658, 0.005296271303249143
                0, 0 
                ];            
                RjcT1 = 0.15 ;
                RjcT2 = 0.16 ;
                RjcD1 = 0.25 ;
                RjcD2 = 0.29 ;
                RjcD5 = 0.31 ;
                RcsT1 = 0.031;      
                RcsT2 = 0.035;      
                RcsD1 = 0.043;      
                RcsD2 = 0.039;  
                RcsD5 = 0.049;  
                weight = 0.398 ;

      end
      
    Semicond.Param.eon_t1=sort(eon_t1);
    Semicond.Param.eon_t2=sort(eon_t2);
    Semicond.Param.eoff_t1=sort(eoff_t1);
    Semicond.Param.eoff_t2=sort(eoff_t2);
    
    Semicond.Param.err_d5=sort(err_d5);
    Semicond.Param.err_d1=sort(err_d1);

    eon_t1x= pchip(eon_t1(:,1),eon_t1(:,2),0:600);
    eon_t2x= pchip(eon_t2(:,1),eon_t2(:,2),0:600);
    
    eoff_t1x= pchip(eoff_t1(:,1),eoff_t1(:,2),0:600);
    eoff_t2x= pchip(eoff_t2(:,1),eoff_t2(:,2),0:600);
    
    err_d5x= pchip(err_d5(:,1),err_d5(:,2),0:600);
    err_d1x= pchip(err_d1(:,1),err_d1(:,2),0:600);

    Semicond.Param.eon_t1x=eon_t1x;
    Semicond.Param.eon_t2x=eon_t2x;
    
    Semicond.Param.eoff_t1x=eoff_t1x;
    Semicond.Param.eoff_t2x=eoff_t2x;
    
    Semicond.Param.err_d5x=err_d5x;
    Semicond.Param.err_d1x=err_d1x;
    
    Semicond.Param.RjcT1 = RjcT1;
    Semicond.Param.RjcT2 = RjcT2;
    Semicond.Param.RjcD1 = RjcD1;
    Semicond.Param.RjcD2 = RjcD2;
    Semicond.Param.RjcD5 = RjcD5;
    Semicond.Param.RcsT1 = RcsT1;
    Semicond.Param.RcsT2 = RcsT2;      
    Semicond.Param.RcsD1 = RcsD1;      
    Semicond.Param.RcsD2 = RcsD2;    
    Semicond.Param.RcsD5 = RcsD5;    
    Semicond.Param.weight = weight ;
    
    case '3L_TTYPE'
      switch Type
          
            case 'SEMiX305TMLI12E4B'
                %% Diode 1  and Switch 2
                % IGBT parameters
                Semicond.Param.U__Bt2 = 300;
                Semicond.Param.rt2 = 0.0041;
                Semicond.Param.Vt2 = 0.7;

                % Diode parameters
                Semicond.Param.U__Bd1 = 300;
                Semicond.Param.rd1 = 0.0028;
                Semicond.Param.Vd1 = 1.4;

                % E reverse recovery
                err_d1=[103.7037037037037, 0.0073947368421052664
                186.41975308641975, 0.007789473684210527
                281.48148148148147, 0.00810526315789474
                353.08641975308643, 0.008342105263157894
                458.02469135802465, 0.008421052631578954
                570.3703703703704, 0.008421052631578954
                658.0246913580247, 0.008578947368421054
                698.7654320987654, 0.008421052631578954
                0, 0
                ];

                % E turn off
                eoff_t2=[101.23456790123453, 0.009605263157894735
                162.962962962963, 0.011894736842105263
                237.037037037037, 0.014500000000000006
                318.51851851851853, 0.01813157894736843
                411.11111111111114, 0.023026315789473686
                483.95061728395063, 0.026500000000000006
                564.1975308641975, 0.029342105263157895
                648.148148148148, 0.03202631578947369
                696.2962962962961, 0.03352631578947368
                0, 0
                ];


                % E turn on
                eon_t2=[103.7037037037037, 0.0017105263157894748
                164.19753086419752, 0.0018684210526315817
                243.20987654320987, 0.0020263157894736886
                322.2222222222221, 0.0025000000000000022
                407.40740740740733, 0.0030526315789473693
                480.2469135802468, 0.0036052631578947364
                574.0740740740741, 0.0041578947368421035
                656.7901234567901, 0.004631578947368424
                696.2962962962961, 0.004868421052631584
                0, 0 
                ];
            
                %% Diode 2  and Switch 1
                % IGBT parameters
                Semicond.Param.U__Bt1 = 300;
                Semicond.Param.rt1 = 0.0056;
                Semicond.Param.Vt1 = 0.8;

                % Diode parameters
                Semicond.Param.U__Bd2 = 300;
                Semicond.Param.rd2 = 0.0016;
                Semicond.Param.Vd2 = 1.15;

                % E reverse recovery
                err_d2=[0.3076923076923208, 0
                102.46153846153845, 0.006888888888888889
                180, 0.007555555555555565
                279.6923076923076, 0.007888888888888904
                364.6153846153845, 0.007777777777777786
                476.6153846153845, 0.008222222222222228
                555.3846153846152, 0.008222222222222228
                625.5384615384614, 0.00811111111111111
                696.9230769230768, 0.00811111111111111
                ];

                % E turn off
                eoff_t1=[1.5384615384615472, 0
                103.69230769230768, 0.01122222222222223
                176.30769230769226, 0.014888888888888889
                253.8461538461538, 0.018888888888888893
                330.1538461538461, 0.022333333333333344
                428.6153846153845, 0.02655555555555556
                506.1538461538461, 0.030888888888888896
                572.6153846153844, 0.0348888888888889
                635.3846153846152, 0.03955555555555556
                696.9230769230768, 0.044111111111111115
                ];

                % E turn on
                eon_t1=[103.69230769230768, 0.002555555555555561
                162.76923076923072, 0.002888888888888899
                0.3076923076923208, -0.0001111111111111035
                225.5384615384615, 0.003555555555555562
                316.6153846153845, 0.0050000000000000044
                407.6923076923076, 0.007444444444444448
                465.53846153846143, 0.009222222222222229
                528.3076923076922, 0.010666666666666672
                598.4615384615383, 0.012444444444444452
                647.6923076923076, 0.013666666666666674
                696.9230769230768, 0.014666666666666675
                ];  
            
                RjcT1 = 0.09 ;
                RjcT2 = 0.17 ;
                RjcD1 = 0.16 ;
                RjcD2 = 0.22 ;
                RcsT1 = 0.021;      
                RcsT2 = 0.059;      
                RcsD1 = 0.041;      
                RcsD2 = 0.059;    
                weight = 0.398 ; 

      end
      
    Semicond.Param.eon_t1=sort(eon_t1);
    Semicond.Param.eon_t2=sort(eon_t2);
    Semicond.Param.eoff_t1=sort(eoff_t1);
    Semicond.Param.eoff_t2=sort(eoff_t2);
    
    Semicond.Param.err_d2=sort(err_d2);
    Semicond.Param.err_d1=sort(err_d1);

    eon_t1x= pchip(eon_t1(:,1),eon_t1(:,2),0:600);
    eon_t2x= pchip(eon_t2(:,1),eon_t2(:,2),0:600);
    
    eoff_t1x= pchip(eoff_t1(:,1),eoff_t1(:,2),0:600);
    eoff_t2x= pchip(eoff_t2(:,1),eoff_t2(:,2),0:600);
    
    err_d2x= pchip(err_d2(:,1),err_d2(:,2),0:600);
    err_d1x= pchip(err_d1(:,1),err_d1(:,2),0:600);

    Semicond.Param.eon_t1x=eon_t1x;
    Semicond.Param.eon_t2x=eon_t2x;
    
    Semicond.Param.eoff_t1x=eoff_t1x;
    Semicond.Param.eoff_t2x=eoff_t2x;
    
    Semicond.Param.err_d2x=err_d2x;
    Semicond.Param.err_d1x=err_d1x;
    
    Semicond.Param.RjcT1 = RjcT1;
    Semicond.Param.RjcT2 = RjcT2;
    Semicond.Param.RjcD1 = RjcD1;
    Semicond.Param.RjcD2 = RjcD2;
    Semicond.Param.RcsT1 = RcsT1;
    Semicond.Param.RcsT2 = RcsT2;      
    Semicond.Param.RcsD1 = RcsD1;      
    Semicond.Param.RcsD2 = RcsD2;      
    Semicond.Param.weight = weight ;

end
    
end