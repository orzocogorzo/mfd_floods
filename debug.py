from mfdfloods import MFD

model = MFD(
  "data/vielha_dtm.tif",
  "data/vielha_mannings.tif",
  5,
  1000,
  False
)

riskies = [
  ((320374.056688216, 4730173.02458468), 28, 14, 9),
  ((320365.742923193, 4730178.58168025), 28, 14, 3),
  ((320357.429158171, 4730184.13877582), 28, 14, 3),
  ((320349.115393148, 4730189.69587139), 28, 14, 9),
  ((320340.801628125, 4730195.25296696), 28, 14, 3),
  ((320332.487863103, 4730200.81006252), 28, 14, 9),
  ((320324.17409808, 4730206.36715809), 28, 14, 9),
  ((320315.860333058, 4730211.92425366), 28, 14, 6),
  ((320307.546568035, 4730217.48134923), 28, 14, 9),
  ((320299.232803012, 4730223.03844479), 28, 14, 9),
  ((320290.91903799, 4730228.59554036), 28, 14, 3),
  ((320282.605272967, 4730234.15263593), 28, 14, 9),
  ((320274.291507945, 4730239.7097315), 28, 14, 3),
  ((320265.977742922, 4730245.26682707), 28, 14, 3),
  ((320257.663977899, 4730250.82392263), 28, 14, 9),
  ((320249.350212877, 4730256.3810182), 28, 14, 3),
  ((320241.036447854, 4730261.93811377), 28, 14, 3),
  ((320232.966123304, 4730267.79013584), 28, 14, 9),
  ((320225.93401366, 4730274.89994836), 28, 14, 3),
  ((320218.901904017, 4730282.00976087), 28, 14, 9),
  ((320211.869794374, 4730289.11957338), 28, 14, 9),
  ((320204.83768473, 4730296.22938589), 28, 14, 9),
  ((320190.773465443, 4730310.44901092), 28, 14, 3),
  ((320183.7413558, 4730317.55882343), 28, 14, 9),
  ((320176.709246156, 4730324.66863594), 28, 14, 3),
  ((320169.677136513, 4730331.77844846), 28, 14, 3),
  ((320162.64502687, 4730338.88826097), 28, 14, 9),
  ((320155.612917226, 4730345.99807348), 28, 14, 3),
  ((320148.580807583, 4730353.10788599), 28, 14, 9),
  ((320142.385604774, 4730360.89996223), 28, 14, 3),
  ((320136.838602811, 4730369.22046517), 28, 14, 9),
  ((320131.291600849, 4730377.54096812), 28, 14, 9),
  ((320125.744598887, 4730385.86147106), 28, 14, 3),
  ((320120.197596925, 4730394.181974), 28, 14, 9),
  ((320114.650594963, 4730402.50247695), 28, 14, 9),
  ((320109.103593, 4730410.82297989), 28, 14, 7),
  ((320081.368583189, 4730452.42549461), 28, 14, 3),
  ((320075.821581227, 4730460.74599755), 28, 14, 9),
  ((320070.274579265, 4730469.06650049), 28, 14, 9),
  ((320064.727577302, 4730477.38700344), 28, 14, 3),
  ((320059.18057534, 4730485.70750638), 28, 14, 9),
  ((320053.633573378, 4730494.02800933), 28, 14, 3),
  ((320049.423129567, 4730503.05860338), 28, 14, 3),
  ((320045.633112013, 4730512.3125629), 28, 14, 3),
]

for risky in riskies:
  floods, drafts, speeds = model.drainpaths(*risky)
