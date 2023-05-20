import React from 'react';
import Login from './pages/Login';
import { Navigate, Route, Routes } from "react-router-dom";


const AppRouter = () => {
  return (
    <Routes>
      <Route path="/login" element={<Login/>}></Route>
      <Route path="/*" element={<Navigate to="/login" />} />
    </Routes>
  );
};

export default AppRouter;