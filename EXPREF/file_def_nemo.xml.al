<?xml version="1.0"?>
    <!-- 
============================================================================================================
=                                           output files definition                                        =
=                                            Define your own files                                         =
=                                         put the variables you want...                                    =
============================================================================================================
    -->

   <file_definition type="one_file" name="@expname@_@freq@_@startdate@_@enddate@" sync_freq="10d" min_digits="4">
    
      <file_group id="1ts" output_freq="1ts"  output_level="10" enabled=".TRUE."/> <!-- 1 time step files -->
      <file_group id="1h" output_freq="1h"  output_level="10" enabled=".TRUE."/> <!-- 1h files -->
      <file_group id="2h" output_freq="2h"  output_level="10" enabled=".TRUE."/> <!-- 2h files -->
      <file_group id="3h" output_freq="3h"  output_level="10" enabled=".TRUE."/> <!-- 3h files -->     
      <file_group id="4h" output_freq="4h"  output_level="10" enabled=".TRUE."/> <!-- 4h files -->
      <file_group id="6h" output_freq="6h"  output_level="10" enabled=".TRUE."/> <!-- 6h files -->
     
      <file_group id="1d" output_freq="1d"  output_level="10" enabled=".FALSE."> <!-- 1d files -->
        <file id="file1" name_suffix="_bioscalar" description="pisces sms variables" >
           <field field_ref="tdenit"   name="tdenit"    unit="TgN/yr" operation="instant" > tdenit * 14. * 86400. * 365. / 1e12 </field>
           <field field_ref="tnfix"    name="tnfix"     unit="TgN/yr" operation="instant" > tnfix * 14. * 86400. * 365. / 1e12 </field>
           <field field_ref="tcflx"    name="tcflx"     unit="PgC/yr" operation="instant" > tcflx * -1. * 12. * 86400. * 365. / 1e15 </field>
           <field field_ref="tcflxcum" name="tcflxcum"  unit="PgC"    operation="instant" > tcflxcum * -1. * 12. / 1e15 </field>
           <field field_ref="tcexp"    name="tcexp"     unit="PgC/yr" operation="instant" > tcexp * 12. * 86400. * 365. / 1e15 </field>
           <field field_ref="tintpp"   name="tintpp"    unit="PgC/yr" operation="instant" > tintpp * 12. * 86400. * 365. / 1e15 </field>
           <field field_ref="pno3tot"  name="pno3tot"   unit="umolN"  > pno3tot * 16. / 122. * 1e6 </field>
           <field field_ref="ppo4tot"  name="ppo4tot"   unit="umolP"  > ppo4tot * 1. / 122. * 1e6 </field>
           <field field_ref="psiltot"  name="psiltot"   unit="umolC"  > psiltot * 1e6  </field>
           <field field_ref="palktot"  name="palktot"   unit="umolC"  > palktot * 1e6  </field>
           <field field_ref="pfertot"  name="pfertot"   unit="nmolFe" > pfertot * 1e9  </field>
        </file>
      </file_group>

      <file_group id="3d" output_freq="3d"  output_level="10" enabled=".TRUE."/> <!-- 3d files -->    
      <file_group id="5d" output_freq="5d"  output_level="10" enabled=".TRUE."/>  <!-- 5d files -->   

      <file_group id="1m" output_freq="1mo" output_level="10" enabled=".TRUE."> <!-- real monthly files -->

	<file id="file2" name_suffix="_ptrc_T" description="pisces sms variables" >
          <field field_ref="DIC"      />
          <field field_ref="Alkalini" />
          <field field_ref="O2"       />
          <field field_ref="CaCO3"    />
          <field field_ref="PO4"      />
          <field field_ref="POC"      />
          <field field_ref="Si"       />
          <field field_ref="PHY"      />
          <field field_ref="PHYN"      />
          <field field_ref="PHYP"      />
          <field field_ref="DOC"      />
          <field field_ref="DON"      />
          <field field_ref="DOP"      />
          <field field_ref="PHY2"     />
          <field field_ref="DIAN"     />
          <field field_ref="DIAP"     />
          <field field_ref="ZOO"     />
          <field field_ref="ZOO2"     />
          <field field_ref="DSi"      />
          <field field_ref="Fer"      />
          <field field_ref="BFe"      />
          <field field_ref="GOC"      />
          <field field_ref="GON"      />
          <field field_ref="GOP"      />
          <field field_ref="SFe"      />
          <field field_ref="DFe"      />
          <field field_ref="GSi"      />
          <field field_ref="NFe"      />
          <field field_ref="NCHL"     />
          <field field_ref="DCHL"     />
          <field field_ref="NO3"      />
          <field field_ref="NH4"      />
          <field field_ref="PIC"     />
          <field field_ref="PICN"     />
          <field field_ref="PICP"     />
          <field field_ref="PFe"      />
          <field field_ref="PCHL"     />
          <field field_ref="PON"     />
          <field field_ref="POP"     />
          <field field_ref="LGW"     />
	</file>
	
	<file id="file3" name_suffix="_diad_T" description="additional pisces diagnostics" >
          <field field_ref="Cflx"     />
          <field field_ref="Dpco2"    />
	</file>

      </file_group>
      <file_group id="2m" output_freq="2mo" output_level="10" enabled=".TRUE."/> <!-- real 2m files -->
      <file_group id="3m" output_freq="3mo" output_level="10" enabled=".TRUE."/> <!-- real 3m files -->
      <file_group id="4m" output_freq="4mo" output_level="10" enabled=".TRUE."/> <!-- real 4m files -->
      <file_group id="6m" output_freq="6mo" output_level="10" enabled=".TRUE."/> <!-- real 6m files -->

      <file_group id="1y"  output_freq="1y" output_level="10" enabled=".TRUE."> <!-- real yearly files -->

	<file id="file4" name_suffix="_ptrc_T" description="pisces sms variables" >
          <field field_ref="DIC"      />
          <field field_ref="Alkalini" />
          <field field_ref="O2"       />
          <field field_ref="CaCO3"    />
          <field field_ref="PO4"      />
          <field field_ref="POC"      />
          <field field_ref="Si"       />
          <field field_ref="PHY"      />
          <field field_ref="PHYN"      />
          <field field_ref="PHYP"      />
          <field field_ref="DOC"      />
          <field field_ref="DON"      />
          <field field_ref="DOP"      />
          <field field_ref="PHY2"     />
          <field field_ref="DIAN"     />
          <field field_ref="DIAP"     />
          <field field_ref="ZOO"     />
          <field field_ref="ZOO2"     />
          <field field_ref="DSi"      />
          <field field_ref="Fer"      />
          <field field_ref="BFe"      />
          <field field_ref="GOC"      />
          <field field_ref="GON"      />
          <field field_ref="GOP"      />
          <field field_ref="SFe"      />
          <field field_ref="DFe"      />
          <field field_ref="GSi"      />
          <field field_ref="NFe"      />
          <field field_ref="NCHL"     />
          <field field_ref="DCHL"     />
          <field field_ref="NO3"      />
          <field field_ref="NH4"      />
          <field field_ref="PIC"     />
          <field field_ref="PICN"     />
          <field field_ref="PICP"     />
          <field field_ref="PFe"      />
          <field field_ref="PCHL"     />
          <field field_ref="PON"     />
          <field field_ref="POP"     />
          <field field_ref="LGW"     />
	</file>

	<file id="file5" name_suffix="_diad_T" description="additional pisces diagnostics" >
          <field field_ref="PH"       />
          <field field_ref="CO3"      />
          <field field_ref="CO3sat"   />
          <field field_ref="PAR"      />
          <field field_ref="PPPHYN"    />
          <field field_ref="PPPHYD"   />
          <field field_ref="PPPHYP"   />
          <field field_ref="PPNEWN"   />
          <field field_ref="PPNEWD"   />
          <field field_ref="PPNEWP"   />
          <field field_ref="PBSi"     />
          <field field_ref="PFeN"     />
          <field field_ref="PFeD"     />
          <field field_ref="PFeP"     />
          <field field_ref="xfracal"  />
          <field field_ref="PCAL"     />
          <field field_ref="DCAL"     />
          <field field_ref="GRAZ1"    />
          <field field_ref="GRAZ2"    />
          <field field_ref="EPC100"   />
          <field field_ref="EPFE100"  />
          <field field_ref="EPSI100"  />
          <field field_ref="EPCAL100" />
          <field field_ref="EXPC"  />
          <field field_ref="EXPFE"  />
          <field field_ref="EXPSI"  />
          <field field_ref="EXPCAL"  />
          <field field_ref="REMIN"    />
          <field field_ref="REMINP"    />
          <field field_ref="REMING"    />
          <field field_ref="Cflx"     />
          <field field_ref="Oflx"     />
          <field field_ref="Kg"       />
          <field field_ref="Dpco2"    />
          <field field_ref="Dpo2"     />
          <field field_ref="Heup"     />
          <field field_ref="Irondep"  />
          <field field_ref="Ironsed"  />
          <field field_ref="Ironice"  />
          <field field_ref="HYDR"  />
          <field field_ref="Nfix"     />
          <field field_ref="MuN"      />
          <field field_ref="MuD"      />
          <field field_ref="MuP"      />
          <field field_ref="LNnut"    />
          <field field_ref="LDnut"    />
          <field field_ref="LPnut"    />
          <field field_ref="LNFe"     />
          <field field_ref="LDFe"     />
          <field field_ref="LPFe"     />
          <field field_ref="LNlight"  />
          <field field_ref="LDlight"  />
          <field field_ref="LPlight"  />
          <field field_ref="SIZEN"  />
          <field field_ref="SIZED"  />
          <field field_ref="SIZEP"  />
          <field field_ref="pdust"    />
          <field field_ref="Fe3"      />
          <field field_ref="FeL1"     />
          <field field_ref="TL1"      />
          <field field_ref="Sdenit"   />
          <field field_ref="Totlig"   />
          <field field_ref="FESCAV"   />
          <field field_ref="FECOLL"   />
          <field field_ref="FEPREC"   />
          <field field_ref="REMINF"    />
          <field field_ref="LGWCOLL"   />
          <field field_ref="LPRODP"   />
          <field field_ref="LDETP"    />
          <field field_ref="LPRODZ"   />
          <field field_ref="LPRODZ2"   />
          <field field_ref="LPRODR"   />
          <field field_ref="LIGREM"   />
          <field field_ref="LIGPR"   />
          <field field_ref="FEZOO"   />
          <field field_ref="FEZOO2"   />
          <field field_ref="BACT" />
          <field field_ref="MLD" />
          <field field_ref="QFePSN"  />
          <field field_ref="QFeNO3N"  />
          <field field_ref="QFeRESN"  />
          <field field_ref="QFePSD"  />
          <field field_ref="QFeNO3D"  />
          <field field_ref="QFeRESD"  />
          <field field_ref="QFePSP"  />
          <field field_ref="QFeNO3P"  />
          <field field_ref="QFeRESP"  />
          <field field_ref="LNN"     />
          <field field_ref="LDN"     />
          <field field_ref="LPN"     />
          <field field_ref="QFEMAXN"     />
          <field field_ref="QFEMAXD"     />
          <field field_ref="QFEMAXP"     />
	</file>
        <file id="file6" name_suffix="_bioscalar" description="pisces sms variables" >
          <field field_ref="tdenit"   name="tdenit"   unit="TgN/yr" > tdenit * 14. * 86400. * 365. / 1e12 </field>
          <field field_ref="tnfix"    name="tnfix"    unit="TgN/yr" > tnfix * 14. * 86400. * 365. / 1e12 </field>
          <field field_ref="tcflx"    name="tcflx"    unit="PgC/yr" > tcflx * -1. * 12. * 86400. * 365. / 1e15 </field>
          <field field_ref="tcflxcum" name="tcflxcum" unit="PgC"    > tcflxcum * -1. * 12. / 1e15 </field>
          <field field_ref="tcexp"    name="tcexp"    unit="PgC/yr" > tcexp * 12. * 86400. * 365. / 1e15 </field>
          <field field_ref="tintpp"   name="tintpp"   unit="PgC/yr" > tintpp * 12. * 86400. * 365. / 1e15 </field>
          <field field_ref="pno3tot"  name="pno3tot"  unit="umolN"  > pno3tot * 16. / 122. * 1e6 </field>
          <field field_ref="ppo4tot"  name="ppo4tot"  unit="umolP"  > ppo4tot * 1. / 122. * 1e6 </field>
          <field field_ref="psiltot"  name="psiltot"  unit="umolC"  > psiltot * 1e6  </field>
          <field field_ref="palktot"  name="palktot"  unit="umolC"  > palktot * 1e6  </field>
          <field field_ref="pfertot"  name="pfertot"  unit="nmolFe" > pfertot * 1e9  </field>
        </file>
      </file_group>

      <file_group id="2y"  output_freq="2y" output_level="10" enabled=".TRUE."/> <!-- real 2y files -->
      <file_group id="5y"  output_freq="5y" output_level="10" enabled=".TRUE."/> <!-- real 5y files -->
      <file_group id="10y" output_freq="10y" output_level="10" enabled=".TRUE."/> <!-- real 10y files -->

   </file_definition>

