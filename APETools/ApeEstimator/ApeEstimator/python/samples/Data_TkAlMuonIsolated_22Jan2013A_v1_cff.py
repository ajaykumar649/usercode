import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/E0B93F3D-AE6F-E211-931F-001F296A52A6.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/74548787-CF6F-E211-9B23-001E0B5FA528.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/D8580CB6-B770-E211-BADB-AC162DACE1B8.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/FC2AA490-7C70-E211-A072-00266CFFBC3C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/FEF88E00-F36F-E211-A80D-001B78CE74FE.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/FC021A6B-9070-E211-BA3D-0017A4771030.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/0E312688-7770-E211-8BB6-00266CFFCCB4.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/F8D16657-6E70-E211-8023-00237DA1B34E.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/0E938275-8F70-E211-A1BB-1CC1DE052030.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/E47BD047-6F70-E211-BA93-1CC1DE041F38.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/2E1A4F42-9F70-E211-86E1-AC162DACE1B8.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/E09871C2-8F70-E211-90F2-1CC1DE048FB0.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/E024283F-E06F-E211-BB1C-00266CFFCA1C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/DEFE1D69-7E70-E211-B5C0-1CC1DE056080.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/E27BB1E4-F36F-E211-A523-1CC1DE0590E8.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/DE066B03-8070-E211-840B-0017A4770034.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/E2B9FD81-F46F-E211-9111-1CC1DE046FB0.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/DCBBACD9-B470-E211-9729-0017A4770034.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/5072FA58-8170-E211-B710-AC162DACB208.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/DA1D4202-9270-E211-837D-1CC1DE051118.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/6058EA7D-8F70-E211-ADF2-00237DA14F86.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/D6209364-8A70-E211-A8B8-0017A4770824.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/323161AF-7870-E211-A558-0017A4770824.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/BE5695BF-7870-E211-BFDB-1CC1DE046FB0.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/B2F04850-6F70-E211-BF8F-1CC1DE04DF70.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/BC5F22BE-9070-E211-988B-0025B3E0228C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/3A5EEA76-7070-E211-90D5-1CC1DE046F18.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/BA7EA258-8170-E211-97D1-00237DA13CA2.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/46D13F8D-8F70-E211-8A43-00266CFFCD00.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/B4C98B81-7070-E211-8E1D-00266CFC3B0C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/90E657D4-8F70-E211-BFC0-0025B3E01FC2.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/7E14B157-9C70-E211-B337-0017A4770038.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/64C413A2-7770-E211-9DCE-0017A4770C24.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/0871D82A-2570-E211-AA64-AC162DABBBA0.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/5CFA10C8-7B70-E211-8DB5-1CC1DE051028.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/12EE435D-2870-E211-8159-1CC1DE041F38.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/5A99D796-8A70-E211-B16D-1CC1DE0500F0.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/12FD5833-BB6F-E211-8D40-0017A4770C08.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/5A4B80B8-9070-E211-B1FE-00237DA10D14.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/18410480-BF6F-E211-AAA9-0017A4770C08.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/545B89BA-7B70-E211-91F2-0017A477083C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/1ABFC1FD-6C70-E211-BFB2-1CC1DE047F98.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/1C31BCBB-F56F-E211-802A-00266CFFBF84.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/20639F17-EF6F-E211-BBEB-1CC1DE0590E8.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/462F9464-7E70-E211-A839-00237DA13C5A.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/02B803D1-B770-E211-8287-00266CFFC940.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/306E29F9-C26F-E211-918E-1CC1DE1D036C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/102CB809-B86F-E211-8688-00266CF65AC4.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/308CFD42-9070-E211-A8BC-00237DA14F86.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/10363F4B-E06F-E211-A3D9-00266CFFBF38.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/08039FC6-6C70-E211-897C-1CC1DE1CEDB2.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/2CDBB71D-9C70-E211-B4F0-0017A4770C1C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/226CAFCE-B770-E211-83A2-00266CFFC550.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/2405922C-9B70-E211-8A69-0017A4770C04.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/4E155B7B-B770-E211-B58F-00266CFFBF64.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/1C5DB126-9B70-E211-B8DF-0017A4771000.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/4E59AFB0-B770-E211-AED0-00266CFEFDE0.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/30000/1696CD0B-8070-E211-8CBF-0017A4770818.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/082B1C96-F36F-E211-AAFD-AC162DABAF78.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/5C0BF9B6-6B70-E211-BEDA-0022649F01AA.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/6494B06F-ED6F-E211-827F-00266CFFC598.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/08342E69-EC6F-E211-B6FA-001F29C4616E.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/FACCB53A-D66F-E211-AEF7-0017A4770034.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/6E36AC79-F06F-E211-962B-1CC1DE0590E8.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/FA960843-CF6F-E211-8849-00237DA13CA2.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/0A3BDABF-B770-E211-92C8-1CC1DE1D16E6.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/F6707F71-F46F-E211-AD1C-00266CFFBF84.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/2011F562-E66F-E211-BC07-0017A4770028.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/F6542291-F36F-E211-B9BC-001B78CE74FE.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/2A0AAFB0-B770-E211-8824-00266CFEFDE0.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/F636A0FF-1A70-E211-B3B0-AC162DACC328.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/2A5E6522-F16F-E211-981F-00266CFFBF80.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/F4FDCB99-B770-E211-9BD7-00266CFFBE68.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/3419B503-1970-E211-AC11-1CC1DE1D014A.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/F4BE9334-AE6F-E211-926C-0017A4770C3C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/78574449-D66F-E211-92A6-00266CFFCD6C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/F4A4B802-2F70-E211-AA2B-001E0B48C1FC.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/80E94B37-F16F-E211-90D9-001F296A4DCC.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/F4720E79-B770-E211-B8C9-00266CFF0034.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/2C7873E5-E06F-E211-A80E-001E0B5FE542.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/F2CEBE10-C66F-E211-8E68-1CC1DE04FF50.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/2CD0230C-C66F-E211-953E-00266CFFC598.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/F2B81356-EC6F-E211-B810-001F296A4DCC.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/2EEC32A6-B770-E211-A843-1CC1DE1CDF2A.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/F21A8B97-B770-E211-AD43-00266CFF0ACC.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/30205596-1A70-E211-BA71-AC162DACC328.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/F007E46B-B870-E211-834F-1CC1DE1D036C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/32636F1E-1A70-E211-93EE-00266CFFC6CC.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/EC203EEF-1670-E211-BD91-78E7D159B710.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/840B1E6E-B770-E211-B427-001CC4164B40.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/EAA64BAC-B96F-E211-BBFA-00266CFF0138.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/3290B5C0-CD6F-E211-8AEB-00237DA13CA2.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/E84FCA3D-2570-E211-BFB3-1CC1DE1D03EA.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/32E53272-1C70-E211-8310-1CC1DE1D014A.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/E6963AC3-6C70-E211-96F6-1CC1DE1CEAEE.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/86FCE7F8-2E70-E211-94B3-00266CFFCB7C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/8A0CECE4-F36F-E211-B10B-00266CFFBF84.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/3E10030C-6C70-E211-BFF3-1CC1DE1D19FE.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/E2625E5A-6C70-E211-A1F3-1CC1DE1CEAEE.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/8ABE038E-BD6F-E211-A2AE-0022649F01AA.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/3EB02938-E36F-E211-8033-00266CFFCA1C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/8E06EEFD-B76F-E211-91D7-00266CFFCD14.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/9C953AE7-D06F-E211-90BB-001E0B5FA528.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/D6F01588-F16F-E211-B253-00237DA16C5E.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/409B1ADE-2270-E211-A95C-1CC1DE1CED1C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/D65111C3-CD6F-E211-85B9-00266CFFBCFC.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/A0EDDA5C-E66F-E211-B537-00237DA46FF4.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/D643946D-1C70-E211-B0C2-00266CFFCD60.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/B0EFA9EF-B770-E211-BCE7-001CC4C10E02.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/D0ED8FEA-2270-E211-9CF1-00266CFFC550.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/B613ADF5-D06F-E211-B07D-1CC1DE1CED1C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/CC12F6AE-F56F-E211-86D1-00266CFFC664.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/40B593C4-B770-E211-8D5A-1CC1DE1CEDB2.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/C8CBFB2E-E36F-E211-87AA-001E0B5FA528.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/BA501B16-1A70-E211-A294-00266CFFCD60.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/C8A6E389-1A70-E211-B867-00266CFFCD60.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/C2E94478-ED6F-E211-9E66-0022649C40A8.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/C4CECB76-BD6F-E211-8070-0017A4770C08.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/44A2754F-C46F-E211-8305-0017A4770830.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/BC60144A-AE6F-E211-87AA-00266CFFC940.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/44B23D7C-2870-E211-A34A-1CC1DE047FD8.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/48790421-F26F-E211-A076-00266CFFBF84.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/4C7E3F58-6C70-E211-97CA-0022649F01AA.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/4E17AF57-2670-E211-B25C-1CC1DE056008.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/9E3526FD-1870-E211-9274-1CC1DE1D023A.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/4EC65668-F16F-E211-B8DB-AC162DABAF78.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/9CEDB716-1770-E211-B37C-1CC1DE1D014A.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/5ABA9931-2670-E211-BBAD-1CC1DE1D03EA.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/92CA9C45-AE6F-E211-83E7-78E7D159B710.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/5C459CDD-E06F-E211-BD43-00266CFFC76C.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/90183BB9-C36F-E211-A209-00266CFFC544.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/8AF5DAF9-1A70-E211-87C3-00266CFFCD60.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/608BFEEF-C26F-E211-B7FC-1CC1DE054078.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/8877D6A8-B96F-E211-BEE8-0017A4771004.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/86F5F873-F06F-E211-85EB-AC162DABAF78.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/6202B117-F26F-E211-9BE7-1CC1DE0590E8.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/680B48C2-BF6F-E211-9E9F-0022649F01AA.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/7C88D676-B96F-E211-B32C-00266CFFBCFC.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/68D5E562-EE6F-E211-8982-AC162DABAF78.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/7003FFFF-F26F-E211-A92B-00266CFFC598.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/6CC01BF2-1670-E211-9C66-1CC1DE1D1FE6.root',
       '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/6EACBF39-CF6F-E211-8D1B-00266CFFBCFC.root' ] );


secFiles.extend( [
               ] )


