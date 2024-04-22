import "ketcher-react/dist/index.css";

import "./App.css";
import ChemEditor from "./pages/chemeditor/ChemEditor";
import React from "react";
import AppNavbar from "./components/AppNavbar";

function App() {
  return (
    <div className="app">
      <AppNavbar />
      <ChemEditor />
    </div>
  );
}

export default App;
