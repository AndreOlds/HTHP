package Fluids
extends Modelica.Icons.Package;

//Annotation
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Line(points = {{-76, -80}, {-62, -30}, {-32, 40}, {4, 66}, {48, 66}, {73, 45}, {62, -8}, {48, -50}, {38, -80}}, color = {64, 64, 64}, smooth = Smooth.Bezier), Line(points = {{-40, 20}, {68, 20}}, color = {175, 175, 175}), Line(points = {{-40, 20}, {-44, 88}, {-44, 88}}, color = {175, 175, 175}), Line(points = {{68, 20}, {86, -58}}, color = {175, 175, 175}), Line(points = {{-60, -28}, {56, -28}}, color = {175, 175, 175}), Line(points = {{-60, -28}, {-74, 84}, {-74, 84}}, color = {175, 175, 175}), Line(points = {{56, -28}, {70, -80}}, color = {175, 175, 175}), Line(points = {{-76, -80}, {38, -80}}, color = {175, 175, 175}), Line(points = {{-76, -80}, {-94, -16}, {-94, -16}}, color = {175, 175, 175})}));


  package Water
    extends ExternalMedia.Media.CoolPropMedium(mediumName = "Water", substanceNames = {"Water"});
  end Water;

  package R1233zdE
    extends ExternalMedia.Media.CoolPropMedium(mediumName = "R1233zdE", substanceNames = {"R1233zdE"});
  end R1233zdE;

  package R1234zeE
    extends ExternalMedia.Media.CoolPropMedium(mediumName = "R1234zeE", substanceNames = {"R1234zeE"});
  end R1234zeE;

  package R601a
    extends ExternalMedia.Media.CoolPropMedium(mediumName = "R601a", substanceNames = {"R601a"});
  end R601a;
  
  package R601
    extends ExternalMedia.Media.CoolPropMedium(mediumName = "R601", substanceNames = {"R601"});
  end R601;
end Fluids;
