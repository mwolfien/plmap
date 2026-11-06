"""Bundle a COMBINE/OMEX archive: SBML-qual + MaBoSS .bnd/.cfg + SED-ML template + metadata.
The .omex is the publication-grade FAIR artifact (one file, self-describing via its manifest)."""
import zipfile, os, datetime

today = datetime.date.today().isoformat()

# ---- SED-ML (L1V4) template. MaBoSS .cfg is the authoritative sim setup; this mirrors its time course.
# MaBoSS config: max_time=40, time_tick=0.5 -> outputEndTime=40, numberOfSteps=80. sample_count (stochastic
# trajectories) is a MaBoSS-native setting with no SED-ML equivalent; it stays in the .cfg.
sedml = '''<?xml version="1.0" encoding="UTF-8"?>
<!-- SED-ML template for the Prevention & Longevity Map coupled model.
     NOTE: MaBoSS is executed from MUS_P1_coupled.bnd/.cfg (the authoritative simulation setup, bundled
     in this archive). Running this SED-ML against MaBoSS goes through BioSimulators/the MaBoSS engine.
     The algorithm KISAO id below (BKMC, Boolean kinetic Monte-Carlo) should be confirmed against the
     current BioSimulators MaBoSS entry before submission. -->
<sedML xmlns="http://sed-ml.org/sed-ml/level1/version4" level="1" version="4">
  <listOfModels>
    <model id="PLM" language="urn:sedml:language:sbml.level-3.version-1" source="PLM_model.sbml"/>
  </listOfModels>
  <listOfSimulations>
    <uniformTimeCourse id="sim" initialTime="0" outputStartTime="0" outputEndTime="40" numberOfSteps="80">
      <algorithm kisaoID="KISAO:0000581"/> <!-- BKMC (MaBoSS); verify term -->
    </uniformTimeCourse>
  </listOfSimulations>
  <listOfTasks>
    <task id="task_main" modelReference="PLM" simulationReference="sim"/>
  </listOfTasks>
  <listOfDataGenerators>
    <dataGenerator id="dg_time"><listOfVariables>
      <variable id="t" symbol="urn:sedml:symbol:time" taskReference="task_main"/></listOfVariables>
      <math xmlns="http://www.w3.org/1998/Math/MathML"><ci>t</ci></math></dataGenerator>
    <dataGenerator id="dg_aging"><listOfVariables>
      <variable id="v_aging" taskReference="task_main"
        target="/sbml:sbml/sbml:model/qual:listOfQualitativeSpecies/qual:qualitativeSpecies[@qual:id='Aging_index']"/>
      </listOfVariables><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>v_aging</ci></math></dataGenerator>
    <dataGenerator id="dg_sarco"><listOfVariables>
      <variable id="v_sarco" taskReference="task_main"
        target="/sbml:sbml/sbml:model/qual:listOfQualitativeSpecies/qual:qualitativeSpecies[@qual:id='Sarcopenia']"/>
      </listOfVariables><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>v_sarco</ci></math></dataGenerator>
  </listOfDataGenerators>
  <listOfOutputs>
    <report id="report_main">
      <listOfDataSets>
        <dataSet id="ds_t" label="time" dataReference="dg_time"/>
        <dataSet id="ds_aging" label="Aging_index" dataReference="dg_aging"/>
        <dataSet id="ds_sarco" label="Sarcopenia" dataReference="dg_sarco"/>
      </listOfDataSets>
    </report>
  </listOfOutputs>
</sedML>
'''
open('simulation.sedml', 'w').write(sedml)

manifest = '''<?xml version="1.0" encoding="UTF-8"?>
<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">
  <content location="." format="http://identifiers.org/combine.specifications/omex"/>
  <content location="./PLM_model.sbml" format="http://identifiers.org/combine.specifications/sbml.level-3.version-1.qual.version-1" master="true"/>
  <content location="./simulation.sedml" format="http://identifiers.org/combine.specifications/sed-ml.level-1.version-4"/>
  <content location="./MUS_P1_coupled.bnd" format="http://purl.org/NET/mediatypes/text/plain"/>
  <content location="./MUS_P1_coupled.cfg" format="http://purl.org/NET/mediatypes/text/plain"/>
  <content location="./MUS_module_v3.json" format="http://purl.org/NET/mediatypes/application/json"/>
  <content location="./README.txt" format="http://purl.org/NET/mediatypes/text/plain"/>
  <content location="./metadata.rdf" format="http://identifiers.org/combine.specifications/omex-metadata"/>
</omexManifest>
'''
open('manifest.xml', 'w').write(manifest)

metadata = f'''<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dcterms="http://purl.org/dc/terms/"
         xmlns:vCard="http://www.w3.org/2006/vcard/ns#">
  <rdf:Description rdf:about=".">
    <dcterms:title>Prevention &amp; Longevity Map — coupled P1+MUS Boolean model</dcterms:title>
    <dcterms:description>Mechanistic aging/prevention Boolean model (nutrient-sensing/NAD core + muscle module)
      with a hallmark-of-aging sub-score layer. SBML-qual = signed influence graph + canonical sign-logic;
      MaBoSS .bnd/.cfg = exact weighted-rate executable dynamics; SED-ML = simulation template.</dcterms:description>
    <dcterms:creator><rdf:Description><vCard:fn><vCard:text>Markus Wolfien, Data-Science ZMI</vCard:text></vCard:fn></rdf:Description></dcterms:creator>
    <dcterms:created><dcterms:W3CDTF>{today}</dcterms:W3CDTF></dcterms:created>
  </rdf:Description>
</rdf:RDF>
'''
open('metadata.rdf', 'w').write(metadata)

readme = ('COMBINE/OMEX archive — Prevention & Longevity Map coupled model\n'
          '================================================================\n'
          'PLM_model.sbml     SBML L3 qualitative model: 98 nodes, 188 signed regulatory edges,\n'
          '                   canonical sign-consistent Boolean function (STRUCTURE + sign).\n'
          'MUS_P1_coupled.bnd MaBoSS network (EXACT executable dynamics: weighted transition rates).\n'
          'MUS_P1_coupled.cfg MaBoSS configuration = the authoritative simulation setup\n'
          '                   (max_time=40, time_tick=0.5, sample_count, initial states).\n'
          'simulation.sedml   SED-ML template (time course mirroring the .cfg). MaBoSS<->SED-ML execution\n'
          '                   runs via BioSimulators; confirm the BKMC KISAO term before submission.\n'
          'MUS_module_v3.json PubMed-curated edge evidence (DOIs).\n'
          'metadata.rdf       OMEX metadata (title, creator, date).\n\n'
          'Why two model files: SBML-qual is the interoperable structural model (GINsim/MINERVA/CaSQ\n'
          'ecosystem, SBGN-AF rendering); the MaBoSS files carry the exact weighted dynamics the SBML-qual\n'
          'standard does not represent. Together they make the model both standard and exactly reproducible.\n')
open('README.txt', 'w').write(readme)

files = ['manifest.xml', 'PLM_model.sbml', 'simulation.sedml', 'MUS_P1_coupled.bnd',
         'MUS_P1_coupled.cfg', 'MUS_module_v3.json', 'metadata.rdf', 'README.txt']
with zipfile.ZipFile('PLM_model.omex', 'w', zipfile.ZIP_DEFLATED) as z:
    for f in files:
        z.write(f)
print('wrote PLM_model.omex with', files)
print('size', os.path.getsize('PLM_model.omex'), 'bytes')
